from copy import deepcopy as _deepcopy

import numpy as _np
import pandas as _pd
import scipy
from scipy import integrate as _integrate

from atmPy.aerosols.size_distribution import moments as _sizedist_moment_conversion
from atmPy.general import timeseries as _timeseries
from atmPy.general import vertical_profile as _vertical_profile
from atmPy.radiation.mie_scattering import bhmie as _bhmie
# import atmPy.aerosols.size_distribution.sizedistribution as _sizedistribution
from  atmPy.aerosols.size_distribution import sizedistribution as _sizedistribution
# import warnings as _warnings

import time
import psutil
import os

# Todo: Docstring is wrong
# todo: This function can be sped up by breaking it apart. Then have OpticalProperties
#       have properties that call the subfunction on demand
def _perform_tmatrixcalculations(ds, wl, n, axis_ratio):
    import pytmatrix.tmatrix
    import pytmatrix.scatter
    rs = ds/2 * 1e-3 #radius in um 

    scattcross = _np.zeros(ds.shape[0])
    for e,r in enumerate(rs):
        tparticle = pytmatrix.tmatrix.Scatterer(radius=r, wavelength= wl*1e-3, m=n, axis_ratio=axis_ratio)
        #### FIXME below is not correct, what you need is a random orientation
        # of particles, that is of coarse more than just paralle and perpendicular,
        # but one mor orientation ... there is an option to integrate over all 
        # orientations, see manual
        scattcross[e] = (pytmatrix.scatter.sca_xsect(tparticle, h_pol = True) + pytmatrix.scatter.sca_xsect(tparticle, h_pol = False))/2 #this will get the average for both polarizations
        print('.', end = '')

    out = _pd.DataFrame({'scattering_crossection': scattcross}, index = rs*2)
    return out

def size_dist2optical_properties(op, sd, aod=False, noOfAngles=100):
    """
    !!!Tis Docstring need fixn
    Calculates the extinction crossection, AOD, phase function, and asymmetry Parameter for each layer.
    plotting the layer and diameter dependent extinction coefficient gives you an idea what dominates the overall AOD.

    Parameters
    ----------
    wavelength: float.
        wavelength of the scattered light, unit: nm
    n: float.
        Index of refraction of the scattering particles

    noOfAngles: int, optional.
        Number of scattering angles to be calculated. This mostly effects calculations which depend on the phase
        function.

    Returns
    -------
    OpticalProperty instance

    """

    # if not _np.any(sd.index_of_refraction):
    #     txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #     raise ValueError(txt)
    # if not sd.sup_optical_properties_wavelength:
    #     txt = 'Please provied wavelength by setting the attribute sup_optical_properties_wavelength (in nm).'
    #     raise AttributeError(txt)

    sd.parameters4reductions._check_opt_prop_param_exist()
    wavelength = sd.parameters4reductions.wavelength.value
    n = sd.parameters4reductions.refractive_index.value
    mie_result = sd.parameters4reductions.mie_result.value
    asphericity = sd.parameters4reductions.asphericity.value
    
    out = {}
    sdls = sd.convert2numberconcentration()
    index = sdls.data.index
    dist_class = type(sdls).__name__

    if dist_class not in ['SizeDist','SizeDist_TS','SizeDist_LS']:
        raise TypeError('this distribution class (%s) can not be converted into optical property yet!'%dist_class)

    # determin if refractive index is changing  or if it is constant
    if isinstance(n, _pd.DataFrame):
        n_multi = True
    else:
        n_multi = False
    if not n_multi:
        if isinstance(mie_result, type(None)):
            # print('do mie', end=' ')
            if asphericity ==1:
                mie, angular_scatt_func = _perform_Miecalculations(_np.array(sdls.bincenters / 1000.), wavelength / 1000., n,
                                                               noOfAngles=noOfAngles)    
                out['angular_scatt_func'] = angular_scatt_func
                mie_result = {'mie': mie,
                              'angular_scatt_func': angular_scatt_func}
            else:
                tmatrix = _perform_tmatrixcalculations(sdls.bincenters, wavelength, n, asphericity)
                mie = tmatrix
        else:
            # print('no need for mie', end=' ')
            if asphericity != 1:
                assert(False), 'currently only a constant refractive index can be considered for non-spherical particles'
            mie = mie_result['mie']
            angular_scatt_func = mie_result['angular_scatt_func']


        # if asphericity == 1:
            
            # mie['angular_scatt_func'] = angular_scatt_func
            
        out['mie_result'] = mie_result #{'mie': mie}
    else:
        out['mie_result'] = None #when n changes there is no single mie_result!
        
    if aod:
        #todo: use function that does a the interpolation instead of the sum?!? I guess this can lead to errors when layers are very thick, since centers are used instea dof edges?
        AOD_layer = _np.zeros((len(sdls.layercenters)))
        
    # The arrays below used to be float32. This was very dangerous, as float64 is the standard and some libaries 
    # as scipy.optimize (the fitting routines) will produce strange results if they see a different float!!!
    extCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)), dtype= _np.float64)
    scattCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)), dtype= _np.float64)
    absCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)), dtype= _np.float64)

    angular_scatt_func_effective = _pd.DataFrame()
    asymmetry_parameter_LS = _np.zeros((len(sdls.data.index.values)))

    #calculate optical properties for each line in the dataFrame
    for i, lc in enumerate(sdls.data.index.values):
        laydata = sdls.data.iloc[i].values # picking a size distribution (either a layer or a point in time)

        if n_multi:
            mie, angular_scatt_func = _perform_Miecalculations(_np.array(sdls.bincenters / 1000.), wavelength / 1000., n.iloc[i].values[0],
                                                               noOfAngles=noOfAngles)
        scattering_coefficient = _get_coefficients(mie.scattering_crossection, laydata)
        scattCoeffPerLayer[i] = scattering_coefficient
        # if asphericity == 1:
        if 'extinction_crossection' in mie:
            extinction_coefficient = _get_coefficients(mie.extinction_crossection, laydata)
            extCoeffPerLayer[i] = extinction_coefficient
        if 'absorption_crossection' in mie:
            absorption_coefficient = _get_coefficients(mie.absorption_crossection, laydata)
            absCoeffPerLayer[i] = absorption_coefficient
        # out['test.extcross'] = mie.extinction_crossection.copy()
        # out['test.extcoeff'] = extinction_coefficient.copy()
        # out['test.laydata'] = laydata

            if aod:
                layerThickness = sdls.layerbounderies[i][1] - sdls.layerbounderies[i][0]
                AOD_perBin = extinction_coefficient * layerThickness
                AOD_layer[i] = AOD_perBin.values.sum()


            scattering_cross_eff = laydata * mie.scattering_crossection
    
            pfe = (laydata * angular_scatt_func).sum(axis=1)  # sum of all angular_scattering_intensities
    
            x_2p = pfe.index.values
            y_2p = pfe.values
    
            # limit to [0,pi]
            y_1p = y_2p[x_2p < _np.pi]
            x_1p = x_2p[x_2p < _np.pi]
    
            y_phase_func = y_1p * 4 * _np.pi / scattering_cross_eff.sum()
            asymmetry_parameter_LS[i] = .5 * _integrate.simps(_np.cos(x_1p) * y_phase_func * _np.sin(x_1p), x_1p)
            angular_scatt_func_effective[
                lc] = pfe * 1e-12 * 1e6  # equivalent to extCoeffPerLayer # similar to  _get_coefficients (converts everthing to meter)
    
        if aod:
            out['AOD'] = AOD_layer[~ _np.isnan(AOD_layer)].sum()
            out['AOD_layer'] = _pd.DataFrame(AOD_layer, index=sdls.layercenters, columns=['AOD per Layer'])
            out['AOD_cum'] = out['AOD_layer'].iloc[::-1].cumsum().iloc[::-1]

    extCoeff_perrow_perbin = _pd.DataFrame(extCoeffPerLayer, index=index, columns=sdls.data.columns)
    scattCoeff_perrow_perbin = _pd.DataFrame(scattCoeffPerLayer, index=index, columns=sdls.data.columns)
    absCoeff_perrow_perbin = _pd.DataFrame(absCoeffPerLayer, index=index, columns=sdls.data.columns)

    # if dist_class == 'SizeDist_TS':
    #     out['extCoeff_perrow_perbin'] = timeseries.TimeSeries_2D(extCoeff_perrow_perbin)
    # if dist_class == 'SizeDist':
    #     out['extCoeff_perrow_perbin'] = _timeseries.TimeSeries(extCoeff_perrow_perbin)
    #     out['scattCoeff_perrow_perbin'] = _timeseries.TimeSeries(scattCoeff_perrow_perbin)
    #     out['absCoeff_perrow_perbin'] = _timeseries.TimeSeries(absCoeff_perrow_perbin)
    # else:
    out['extCoeff_perrow_perbin'] = extCoeff_perrow_perbin
    out['scattCoeff_perrow_perbin'] = scattCoeff_perrow_perbin
    out['absCoeff_perrow_perbin'] = absCoeff_perrow_perbin

    out['parent_type'] = dist_class
    out['asymmetry_param'] = _pd.DataFrame(asymmetry_parameter_LS, index=index,
                                           columns=['asymmetry_param'])

    out['wavelength'] = wavelength
    out['index_of_refraction'] = n
    out['bin_centers'] = sdls.bincenters
    out['bins'] = sdls.bins
    out['binwidth'] = sdls.binwidth
    out['distType'] = sdls.distributionType
    out['angular_scatt_func'] = angular_scatt_func_effective.transpose()

    ### test values
    # out['mie_curve_ext'] = mie.extinction_crossection
    # out['mie_inst'] = mie
    return out
    
def size_dist2optical_properties_old(op, sd, aod=False, noOfAngles=100):
    """
    !!!Tis Docstring need fixn
    Calculates the extinction crossection, AOD, phase function, and asymmetry Parameter for each layer.
    plotting the layer and diameter dependent extinction coefficient gives you an idea what dominates the overall AOD.

    Parameters
    ----------
    wavelength: float.
        wavelength of the scattered light, unit: nm
    n: float.
        Index of refraction of the scattering particles

    noOfAngles: int, optional.
        Number of scattering angles to be calculated. This mostly effects calculations which depend on the phase
        function.

    Returns
    -------
    OpticalProperty instance

    """

    # if not _np.any(sd.index_of_refraction):
    #     txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #     raise ValueError(txt)
    # if not sd.sup_optical_properties_wavelength:
    #     txt = 'Please provied wavelength by setting the attribute sup_optical_properties_wavelength (in nm).'
    #     raise AttributeError(txt)

    sd.parameters4reductions._check_opt_prop_param_exist()
    wavelength = sd.parameters4reductions.wavelength.value
    n = sd.parameters4reductions.refractive_index.value
    mie_result = sd.parameters4reductions.mie_result.value
    out = {}
    sdls = sd.convert2numberconcentration()
    index = sdls.data.index
    dist_class = type(sdls).__name__

    if dist_class not in ['SizeDist','SizeDist_TS','SizeDist_LS']:
        raise TypeError('this distribution class (%s) can not be converted into optical property yet!'%dist_class)

    # determin if refractive index is changing  or if it is constant
    if isinstance(n, _pd.DataFrame):
        n_multi = True
    else:
        n_multi = False
    if not n_multi:
        if isinstance(mie_result, type(None)):
            # print('do mie', end=' ')
            mie, angular_scatt_func = _perform_Miecalculations(_np.array(sdls.bincenters / 1000.), wavelength / 1000., n,
                                                               noOfAngles=noOfAngles)            
        else:
            # print('no need for mie', end=' ')
            # if asphericity != 1:
            #     assert(False), 'currently only a constant refractive index can be considered for non-spherical particles'
            mie = mie_result['mie']
            angular_scatt_func = mie_result['angular_scatt_func']

        out['mie_result'] = {'mie': mie, 'angular_scatt_func': angular_scatt_func}

    if aod:
        #todo: use function that does a the interpolation instead of the sum?!? I guess this can lead to errors when layers are very thick, since centers are used instea dof edges?
        AOD_layer = _np.zeros((len(sdls.layercenters)))

    extCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)), dtype= _np.float32)
    scattCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)), dtype= _np.float32)
    absCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)), dtype= _np.float32)

    angular_scatt_func_effective = _pd.DataFrame()
    asymmetry_parameter_LS = _np.zeros((len(sdls.data.index.values)))

    #calculate optical properties for each line in the dataFrame
    for i, lc in enumerate(sdls.data.index.values):
        laydata = sdls.data.iloc[i].values # picking a size distribution (either a layer or a point in time)

        if n_multi:
            mie, angular_scatt_func = _perform_Miecalculations(_np.array(sdls.bincenters / 1000.), wavelength / 1000., n.iloc[i].values[0],
                                                               noOfAngles=noOfAngles)
        extinction_coefficient = _get_coefficients(mie.extinction_crossection, laydata)
        scattering_coefficient = _get_coefficients(mie.scattering_crossection, laydata)
        absorption_coefficient = _get_coefficients(mie.absorption_crossection, laydata)

        out['test.extcross'] = mie.extinction_crossection.copy()
        out['test.extcoeff'] = extinction_coefficient.copy()
        out['test.laydata'] = laydata

        if aod:
            layerThickness = sdls.layerbounderies[i][1] - sdls.layerbounderies[i][0]
            AOD_perBin = extinction_coefficient * layerThickness
            AOD_layer[i] = AOD_perBin.values.sum()

        extCoeffPerLayer[i] = extinction_coefficient
        scattCoeffPerLayer[i] = scattering_coefficient
        absCoeffPerLayer[i] = absorption_coefficient

        scattering_cross_eff = laydata * mie.scattering_crossection

        pfe = (laydata * angular_scatt_func).sum(axis=1)  # sum of all angular_scattering_intensities

        x_2p = pfe.index.values
        y_2p = pfe.values

        # limit to [0,pi]
        y_1p = y_2p[x_2p < _np.pi]
        x_1p = x_2p[x_2p < _np.pi]

        y_phase_func = y_1p * 4 * _np.pi / scattering_cross_eff.sum()
        asymmetry_parameter_LS[i] = .5 * _integrate.simps(_np.cos(x_1p) * y_phase_func * _np.sin(x_1p), x_1p)
        angular_scatt_func_effective[
            lc] = pfe * 1e-12 * 1e6  # equivalent to extCoeffPerLayer # similar to  _get_coefficients (converts everthing to meter)

    if aod:
        out['AOD'] = AOD_layer[~ _np.isnan(AOD_layer)].sum()
        out['AOD_layer'] = _pd.DataFrame(AOD_layer, index=sdls.layercenters, columns=['AOD per Layer'])
        out['AOD_cum'] = out['AOD_layer'].iloc[::-1].cumsum().iloc[::-1]

    extCoeff_perrow_perbin = _pd.DataFrame(extCoeffPerLayer, index=index, columns=sdls.data.columns)
    scattCoeff_perrow_perbin = _pd.DataFrame(scattCoeffPerLayer, index=index, columns=sdls.data.columns)
    absCoeff_perrow_perbin = _pd.DataFrame(absCoeffPerLayer, index=index, columns=sdls.data.columns)

    # if dist_class == 'SizeDist_TS':
    #     out['extCoeff_perrow_perbin'] = timeseries.TimeSeries_2D(extCoeff_perrow_perbin)
    # if dist_class == 'SizeDist':
    #     out['extCoeff_perrow_perbin'] = _timeseries.TimeSeries(extCoeff_perrow_perbin)
    #     out['scattCoeff_perrow_perbin'] = _timeseries.TimeSeries(scattCoeff_perrow_perbin)
    #     out['absCoeff_perrow_perbin'] = _timeseries.TimeSeries(absCoeff_perrow_perbin)
    # else:
    out['extCoeff_perrow_perbin'] = extCoeff_perrow_perbin
    out['scattCoeff_perrow_perbin'] = scattCoeff_perrow_perbin
    out['absCoeff_perrow_perbin'] = absCoeff_perrow_perbin

    out['parent_type'] = dist_class
    out['asymmetry_param'] = _pd.DataFrame(asymmetry_parameter_LS, index=index,
                                           columns=['asymmetry_param'])

    out['wavelength'] = wavelength
    out['index_of_refraction'] = n
    out['bin_centers'] = sdls.bincenters
    out['bins'] = sdls.bins
    out['binwidth'] = sdls.binwidth
    out['distType'] = sdls.distributionType
    out['angular_scatt_func'] = angular_scatt_func_effective.transpose()

    ### test values
    out['mie_curve_ext'] = mie.extinction_crossection
    out['mie_inst'] = mie
    return out


def DEPRECATED_size_dist2optical_properties(sd, aod=False, noOfAngles=100):
    """
    !!!Tis Docstring need fixn
    Calculates the extinction crossection, AOD, phase function, and asymmetry Parameter for each layer.
    plotting the layer and diameter dependent extinction coefficient gives you an idea what dominates the overall AOD.

    Parameters
    ----------
    wavelength: float.
        wavelength of the scattered light, unit: nm
    n: float.
        Index of refraction of the scattering particles

    noOfAngles: int, optional.
        Number of scattering angles to be calculated. This mostly effects calculations which depend on the phase
        function.

    Returns
    -------
    OpticalProperty instance

    """

    # if not _np.any(sd.index_of_refraction):
    #     txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #     raise ValueError(txt)
    # if not sd.sup_optical_properties_wavelength:
    #     txt = 'Please provied wavelength by setting the attribute sup_optical_properties_wavelength (in nm).'
    #     raise AttributeError(txt)

    sd.optical_properties_settings._check()
    wavelength = sd.optical_properties_settings.wavelength.value
    n = sd.optical_properties_settings.refractive_index.value
    out = {}
    sdls = sd.convert2numberconcentration()
    index = sdls.data.index
    dist_class = type(sdls).__name__

    if dist_class not in ['SizeDist','SizeDist_TS','SizeDist_LS']:
        raise TypeError('this distribution class (%s) can not be converted into optical property yet!'%dist_class)

    # determin if index of refraction changes or if it is constant
    if isinstance(n, _pd.DataFrame):
        n_multi = True
    else:
        n_multi = False
    if not n_multi:
        mie, angular_scatt_func = _perform_Miecalculations(_np.array(sdls.bincenters / 1000.), wavelength / 1000., n,
                                                           noOfAngles=noOfAngles)

    if aod:
        #todo: use function that does a the interpolation instead of the sum?!? I guess this can lead to errors when layers are very thick, since centers are used instea dof edges?
        AOD_layer = _np.zeros((len(sdls.layercenters)))

    extCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)))
    scattCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)))
    absCoeffPerLayer = _np.zeros((len(sdls.data.index.values), len(sdls.bincenters)))

    angular_scatt_func_effective = _pd.DataFrame()
    asymmetry_parameter_LS = _np.zeros((len(sdls.data.index.values)))

    #calculate optical properties for each line in the dataFrame
    for i, lc in enumerate(sdls.data.index.values):
        laydata = sdls.data.iloc[i].values # picking a size distribution (either a layer or a point in time)

        if n_multi:
            mie, angular_scatt_func = _perform_Miecalculations(_np.array(sdls.bincenters / 1000.), wavelength / 1000., n.iloc[i].values[0],
                                                               noOfAngles=noOfAngles)
        extinction_coefficient = _get_coefficients(mie.extinction_crossection, laydata)
        scattering_coefficient = _get_coefficients(mie.scattering_crossection, laydata)
        absorption_coefficient = _get_coefficients(mie.absorption_crossection, laydata)

        if aod:
            layerThickness = sdls.layerbounderies[i][1] - sdls.layerbounderies[i][0]
            AOD_perBin = extinction_coefficient * layerThickness
            AOD_layer[i] = AOD_perBin.values.sum()

        extCoeffPerLayer[i] = extinction_coefficient
        scattCoeffPerLayer[i] = scattering_coefficient
        absCoeffPerLayer[i] = absorption_coefficient

        scattering_cross_eff = laydata * mie.scattering_crossection

        pfe = (laydata * angular_scatt_func).sum(axis=1)  # sum of all angular_scattering_intensities

        x_2p = pfe.index.values
        y_2p = pfe.values

        # limit to [0,pi]
        y_1p = y_2p[x_2p < _np.pi]
        x_1p = x_2p[x_2p < _np.pi]

        y_phase_func = y_1p * 4 * _np.pi / scattering_cross_eff.sum()
        asymmetry_parameter_LS[i] = .5 * _integrate.simps(_np.cos(x_1p) * y_phase_func * _np.sin(x_1p), x_1p)
        angular_scatt_func_effective[
            lc] = pfe * 1e-12 * 1e6  # equivalent to extCoeffPerLayer # similar to  _get_coefficients (converts everthing to meter)

    if aod:
        out['AOD'] = AOD_layer[~ _np.isnan(AOD_layer)].sum()
        out['AOD_layer'] = _pd.DataFrame(AOD_layer, index=sdls.layercenters, columns=['AOD per Layer'])
        out['AOD_cum'] = out['AOD_layer'].iloc[::-1].cumsum().iloc[::-1]

    extCoeff_perrow_perbin = _pd.DataFrame(extCoeffPerLayer, index=index, columns=sdls.data.columns)
    scattCoeff_perrow_perbin = _pd.DataFrame(scattCoeffPerLayer, index=index, columns=sdls.data.columns)
    absCoeff_perrow_perbin = _pd.DataFrame(absCoeffPerLayer, index=index, columns=sdls.data.columns)

    # if dist_class == 'SizeDist_TS':
    #     out['extCoeff_perrow_perbin'] = timeseries.TimeSeries_2D(extCoeff_perrow_perbin)
    if dist_class == 'SizeDist':
        out['extCoeff_perrow_perbin'] = _timeseries.TimeSeries(extCoeff_perrow_perbin)
        out['scattCoeff_perrow_perbin'] = _timeseries.TimeSeries(scattCoeff_perrow_perbin)
        out['absCoeff_perrow_perbin'] = _timeseries.TimeSeries(absCoeff_perrow_perbin)
    else:
        out['extCoeff_perrow_perbin'] = extCoeff_perrow_perbin
        out['scattCoeff_perrow_perbin'] = scattCoeff_perrow_perbin
        out['absCoeff_perrow_perbin'] = absCoeff_perrow_perbin
    # extCoeff_perrow = pd.DataFrame(extCoeff_perrow_perbin.sum(axis=1), columns=['ext_coeff'])
    # if index.dtype == '<M8[ns]':
    #     out['extCoeff_perrow'] = timeseries.TimeSeries(extCoeff_perrow)
    # else:
    #     out['extCoeff_perrow'] = extCoeff_perrow

    out['parent_type'] = dist_class
    out['asymmetry_param'] = _pd.DataFrame(asymmetry_parameter_LS, index=index,
                                           columns=['asymmetry_param'])
    # out['asymmetry_param_alt'] = pd.DataFrame(asymmetry_parameter_LS_alt, index=sdls.layercenters, columns = ['asymmetry_param_alt'])
    # out['OptPropInstance']= OpticalProperties(out, self.bins)
    out['wavelength'] = wavelength
    out['index_of_refraction'] = n
    out['bin_centers'] = sdls.bincenters
    out['bins'] = sdls.bins
    out['binwidth'] = sdls.binwidth
    out['distType'] = sdls.distributionType
    out['angular_scatt_func'] = angular_scatt_func_effective
    # opt_properties = OpticalProperties(out, self.bins)
    # opt_properties.wavelength = wavelength
    # opt_properties.index_of_refractio = n
    # opt_properties.angular_scatt_func = angular_scatt_func_effective  # This is the formaer phase_fct, but since it is the angular scattering intensity, i changed the name
    # opt_properties.parent_dist_LS = self
    if dist_class == 'SizeDist_TS':
        return OpticalProperties_TS(out, parent = sd)
    elif dist_class == 'SizeDist_LS':
        return OpticalProperties_VP(out, parent= sd)
    return out


def hemispheric_backscattering(osf_df):
    """scattering into backwards hemisphere from angulare scattering intensity

    Parameters
    ----------
    osf_df: pandas DataFrame
        This contains the angulare scattering intensity with column names giving the
        angles in radiant


    Returns
    -------
    pandas data frame with the scattering intensities
    """
    import pdb
    # pdb.set_trace()
    def ang_scat_funk2bs(index,ol):
        x = index #_np.deg2rad(index)
        f = ol
        # pdb.set_trace()
        # my phase function goes all the way to two py
        f = f[x < _np.pi]
        x = x[x < _np.pi]
        f_b = f[x >= _np.pi / 2.]
        x_b = x[x >= _np.pi / 2.]
        # pdb.set_trace()
        res_b = 2 * _np.pi * _integrate.simps(f_b * _np.sin(x_b), x_b)
        return res_b

    bs = _np.zeros(osf_df.shape[0])
    index = osf_df.columns
    for i in range(osf_df.shape[0]):
        ol = osf_df.iloc[i,:].values
        bs[i] = ang_scat_funk2bs(index,ol)
    bs = _pd.DataFrame(bs, index = osf_df.index)
    return bs

def hemispheric_forwardscattering(osf_df):
    """scattering into forward hemisphere from angulare scattering intensity

    Parameters
    ----------
    osf_df: pandas DataFrame
        This contains the angulare scattering intensity with column names giving the
        angles in radiant

    Returns
    -------
    pandas data frame with the scattering intensities
    """

    def ang_scat_funk2fs(index,ol):
        x = index #ol.index.values
        f = ol

        # my phase function goes all the way to two py
        f = f[x < _np.pi]
        x = x[x < _np.pi]
        f_f = f[x < _np.pi / 2.]
        x_f = x[x < _np.pi / 2.]

        res_f = 2 * _np.pi * _integrate.simps(f_f * _np.sin(x_f), x_f)
        return res_f

    fs = _np.zeros(osf_df.shape[0])
    index = osf_df.columns
    for i in range(osf_df.shape[0]):
        ol = osf_df.iloc[i,:].values
        fs[i] = ang_scat_funk2fs(index,ol)
    fs = _pd.DataFrame(fs, index = osf_df.index)
    return fs



#Todo: bins are redundand
# Todo: some functions should be switched of
# todo: right now this for layer and time series, not ok
class OpticalProperties(object):
    def __init__(self, parent):
        self._parent_sizedist = parent
        self.parameters = _sizedistribution._Parameters4Reductions_opt_prop(parent)

        # self.asymmetry_param = data['asymmetry_param']

        self._extinction_coeff = None
        self._scattering_coeff = None
        self._absorption_coeff = None
        self._mie_result = None

        self._hemispheric_backscattering = None
        # self._hemispheric_backscattering_ratio = None
        self._hemispheric_forwardscattering = None
        # self._hemispheric_forwardscattering_ratio = None

        self._optical_porperties_pv = None


        self.mean_effective_diameter = None
        self._parent_type = type(parent).__name__
        self.bins = parent.bins
        self.binwidth = parent.binwidth
        self.distributionType = parent.distributionType
        # self._data_period = self.parent_sizedist._data_period





    
    @property
    def scattering_coeff_per_bin(self):
        self._scattering_coeff_per_bin = self._optical_porperties['scattCoeff_perrow_perbin']
        return self._scattering_coeff_per_bin
    
    @property
    def scattering_coeff(self):
        self._scattering_coeff = _pd.DataFrame(self.scattering_coeff_per_bin.sum(axis=1), columns=['scatt_coeff_m^1'])
        return self._scattering_coeff
    
    @property
    def absorption_coeff_per_bin(self):
        assert(self.parameters.asphericity.value==1), 'only scattering works with assymmetric particles so far. Should not be very hard to implement though!!! fix'
        self._absorption_coeff_per_bin = self._optical_porperties['absCoeff_perrow_perbin']
        return self._absorption_coeff_per_bin
    
    @property
    def absorption_coeff(self):
        self._absorption_coeff = _pd.DataFrame(self.absorption_coeff_per_bin.sum(axis=1), columns=['abs_coeff_m^1'])
        return self._absorption_coeff
    
    @property
    def extinction_coeff_per_bin(self):
        assert(self.parameters.asphericity.value==1), 'only scattering works with assymmetric particles so far. Should not be very hard to implement though!!! fix'
        self._extinction_coeff_per_bin = self._optical_porperties['extCoeff_perrow_perbin']
        return self._extinction_coeff_per_bin
    
    @property
    def extinction_coeff(self):
        self._extinction_coeff = _pd.DataFrame(self.extinction_coeff_per_bin.sum(axis=1), columns=['ext_coeff_m^1'])
        return self._extinction_coeff
    
    @property
    def angular_scatt_func(self):
        assert(self.parameters.asphericity.value==1), 'only scattering works with assymmetric particles so far. Should not be very hard to implement though!!! fix'
        self._angular_scatt_func = self._optical_porperties['angular_scatt_func']
        return self._angular_scatt_func
    
    @property
    def _optical_porperties(self):
        # assert(False), 'un_comment all the stuff below or even better: put them into each property'
        if not self._optical_porperties_pv:
            data = size_dist2optical_properties(self, self._parent_sizedist)
            self._optical_porperties_pv = data

            # self._angular_scatt_func = data['angular_scatt_func']

            ####
            # self.parameters.mie_result = data['mie_result']
            self._mie_result = data['mie_result']
        return self._optical_porperties_pv

    @property
    def mie_result(self):
        self._optical_porperties
        return self._mie_result


    @property
    def backscattering(self):
        """
        Scattering at exactly pi (180°).

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        return self.angular_scatt_func.loc[:,_np.pi]

    @property
    def hemispheric_backscattering(self):
        if not _np.any(self._hemispheric_backscattering):
            self._hemispheric_backscattering = hemispheric_backscattering(self.angular_scatt_func)
            self._hemispheric_backscattering_ratio = _pd.DataFrame(
                self._hemispheric_backscattering.iloc[:, 0] / self._scattering_coeff.iloc[:, 0],
                columns=['hem_back_scatt_ratio'])
        return self._hemispheric_backscattering

    @property
    def hemispheric_backscattering_ratio(self):
        self.hemispheric_backscattering
        # if not _np.any(self._hemispheric_backscattering_ratio):
        #     self._hemispheric_backscattering_ratio = _pd.DataFrame(self.hemispheric_backscattering.iloc[:,0] / self._scattering_coeff.iloc[:,0], columns=['hem_beck_scatt_ratio'])
        return self._hemispheric_backscattering_ratio

    @property
    def hemispheric_forwardscattering(self):
        if not _np.any(self._hemispheric_forwardscattering):
            self._hemispheric_forwardscattering = hemispheric_forwardscattering(self.angular_scatt_func)
            self._hemispheric_forwardscattering_ratio = _pd.DataFrame(self._hemispheric_forwardscattering.iloc[:, 0] /  self._scattering_coeff.iloc[:, 0],
                columns=['hem_forward_scatt_ratio'])
        return self._hemispheric_forwardscattering

    @property
    def hemispheric_forwardscattering_ratio(self):
        self.hemispheric_forwardscattering
        # if not _np.any(self._hemispheric_forwardscattering_ratio):
        #     self._hemispheric_forwardscattering_ratio = self.hemispheric_forwardscattering / self.scattering_coeff
        return self._hemispheric_forwardscattering_ratio

    def convert_between_moments(self, moment, verbose = False):
        return _sizedist_moment_conversion.convert(self,moment, verbose = verbose)

    def copy(self):
        return _deepcopy(self)




#Todo: bins are redundand
# Todo: some functions should be switched of
# todo: right now this for layer and time series, not ok
class DEPRECATEDOpticalProperties(object):
    def __init__(self, data, parent = None):
        self.parent_sizedist = parent

        self.data_orig = data
        self.wavelength =  data['wavelength']
        self.index_of_refraction = data['index_of_refraction']
        self.extinction_coeff_per_bin = data['extCoeff_perrow_perbin']
        self.scattering_coeff_per_bin = data['scattCoeff_perrow_perbin']
        self.absorption_coeff_per_bin = data['absCoeff_perrow_perbin']
        self.angular_scatt_func = data['angular_scatt_func']

        # self.asymmetry_param = data['asymmetry_param']

        self.__extinction_coeff_sum_along_d = None
        self.__scattering_coeff_sum_along_d = None
        self.__absorption_coeff_sum_along_d = None


        self.mean_effective_diameter = None
        self._parent_type = data['parent_type']
        self.bins = data['bins']
        self.binwidth = data['binwidth']
        self.distributionType = data['distType']
        # self._data_period = self.parent_sizedist._data_period


    # @property
    # def mean_effective_diameter(self):
    #     if not self.__mean_effective_diameter:


    # # todo: remove
    # @property
    # def extinction_coeff_sum_along_d(self):
    #     _warnings.warn('extinction_coeff_sum_along_d is deprecated and will be removed in future versions. Use extingction_coeff instead')
    #     if not _np.any(self.__extinction_coeff_sum_along_d):
    #         data = self.extinction_coeff_per_bin.data.sum(axis = 1)
    #         df = _pd.DataFrame()
    #         df['ext_coeff_m^1'] = data
    #         if self._parent_type == 'SizeDist_TS':
    #             self.__extinction_coeff_sum_along_d = _timeseries.TimeSeries(df)
    #         elif self._parent_type == 'SizeDist':
    #             self.__extinction_coeff_sum_along_d = df
    #         else:
    #             raise TypeError('not possible for this distribution type')
    #         self.__extinction_coeff_sum_along_d._data_period = self._data_period
    #     return self.__extinction_coeff_sum_along_d
    #
    # # todo: remove
    # @extinction_coeff_sum_along_d.setter
    # def extinction_coeff_sum_along_d(self, data):
    #     self.__extinction_coeff_sum_along_d = data

    @property
    def extinction_coeff(self):
        if not _np.any(self.__extinction_coeff_sum_along_d):
            data = self.extinction_coeff_per_bin.data.sum(axis=1)
            df = _pd.DataFrame()
            df['ext_coeff_m^1'] = data
            if self._parent_type == 'SizeDist_TS':
                self.__extinction_coeff_sum_along_d = _timeseries.TimeSeries(df)
                self.__extinction_coeff_sum_along_d._data_period = self._data_period
            elif self._parent_type == 'SizeDist_LS':
                self.__extinction_coeff_sum_along_d = _vertical_profile.VerticalProfile(df)
            elif self._parent_type == 'SizeDist':
                self.__extinction_coeff_sum_along_d = df
            else:
                raise TypeError('not possible for this distribution type')
        return self.__extinction_coeff_sum_along_d

    @extinction_coeff.setter
    def extinction_coeff(self, data):
        self.__extinction_coeff_sum_along_d = data

    @property
    def scattering_coeff(self):
        if not _np.any(self.__scattering_coeff_sum_along_d):
            data = self.scattering_coeff_per_bin.data.sum(axis=1)
            df = _pd.DataFrame()
            df['scatt_coeff_m^1'] = data
            if self._parent_type == 'SizeDist_TS':
                self.__scattering_coeff_sum_along_d = _timeseries.TimeSeries(df)
            elif self._parent_type == 'SizeDist':
                self.__scattering_coeff_sum_along_d = df
            else:
                raise TypeError('not possible for this distribution type')
            self.__scattering_coeff_sum_along_d._data_period = self._data_period
        return self.__scattering_coeff_sum_along_d

    @scattering_coeff.setter
    def scattering_coeff(self, data):
        self.__scattering_coeff_sum_along_d = data



    @property
    def absorption_coeff(self):
        if not _np.any(self.__absorption_coeff_sum_along_d):
            data = self.absorption_coeff_per_bin.data.sum(axis=1)
            df = _pd.DataFrame()
            df['abs_coeff_m^1'] = data
            if self._parent_type == 'SizeDist_TS':
                self.__absorption_coeff_sum_along_d = _timeseries.TimeSeries(df)
            elif self._parent_type == 'SizeDist':
                self.__absorption_coeff_sum_along_d = df
            else:
                raise TypeError('not possible for this distribution type')
            self.__absorption_coeff_sum_along_d._data_period = self._data_period
        return self.__absorption_coeff_sum_along_d

    @absorption_coeff.setter
    def absorption_coeff(self, data):
        self.__absorption_coeff_sum_along_d = data


    @property
    def hemispheric_backscattering(self):
        if not self.__hemispheric_backscattering:
            self.__hemispheric_backscattering = hemispheric_backscattering(self.angular_scatt_func)
        return self.__hemispheric_backscattering

    @property
    def hemispheric_forwardscattering(self):
        if not self.__hemispheric_forwardscattering:
            self.__hemispheric_forwardscattering = hemispheric_forwardscattering(self.angular_scatt_func)
        return self.__hemispheric_forwardscattering

    @property
    def hemispheric_backscattering_ratio(self):
        if not self.__hemispheric_backscattering_ratio:
            self.__hemispheric_backscattering_ratio = self.hemispheric_backscattering / self.scattering_coeff
        return self.__hemispheric_backscattering_ratio

    @property
    def hemispheric_forwardscattering_ratio(self):
        if not self.hemispheric_forwardscattering_ratio:
            self.__hemispheric_forwardscattering_ratio = self.hemispheric_forwardscattering / self.scattering_coeff
        return self.__hemispheric_forwardscattering_ratio

    def convert_between_moments(self, moment, verbose = False):
        return _sizedist_moment_conversion.convert(self,moment, verbose = verbose)

    def copy(self):
        return _deepcopy(self)



class OpticalProperties_TS(OpticalProperties):
    @property
    def hemispheric_forwardscattering(self):
        return _timeseries.TimeSeries(super().hemispheric_forwardscattering, sampling_period = self._parent_sizedist._data_period)

    @property
    def hemispheric_backscattering(self):
        return _timeseries.TimeSeries(super().hemispheric_backscattering, sampling_period = self._parent_sizedist._data_period)

    @property
    def hemispheric_backscattering_ratio(self):
        return _timeseries.TimeSeries(super().hemispheric_backscattering_ratio, sampling_period = self._parent_sizedist._data_period)

    @property
    def hemispheric_forwardscattering_ratio(self):
        return _timeseries.TimeSeries(super().hemispheric_forwardscattering_ratio, sampling_period = self._parent_sizedist._data_period)

    @property
    def absorption_coeff(self):
        return _timeseries.TimeSeries(super().absorption_coeff, sampling_period = self._parent_sizedist._data_period)

    @property
    def extinction_coeff(self):
        return _timeseries.TimeSeries(super().extinction_coeff, sampling_period = self._parent_sizedist._data_period)

    @property
    def scattering_coeff(self):
        return _timeseries.TimeSeries(super().scattering_coeff, sampling_period = self._parent_sizedist._data_period)

    # def __init__(self, *args, **kwargs):
    #     super().__init__(*args, **kwargs)
    #     self.extinction_coeff_per_bin = _timeseries.TimeSeries_2D(self.extinction_coeff_per_bin)
    #     self.extinction_coeff_per_bin._data_period = self.parent_sizedist._data_period
    #
    #     self.scattering_coeff_per_bin = _timeseries.TimeSeries_2D(self.scattering_coeff_per_bin)
    #     self.scattering_coeff_per_bin._data_period = self.parent_sizedist._data_period
    #
    #     self.absorption_coeff_per_bin = _timeseries.TimeSeries_2D(self.absorption_coeff_per_bin)
    #     self.absorption_coeff_per_bin._data_period = self.parent_sizedist._data_period
    #
    #     self.angular_scatt_func = _timeseries.TimeSeries_2D(self.angular_scatt_func.transpose())
    #     self.angular_scatt_func._data_period = self.parent_sizedist._data_period
    #
    #     self.__hemispheric_forwardscattering = None
    #     self.__hemispheric_backscattering = None
    #     self.__hemispheric_backscattering_ratio = None
    #     self.__hemispheric_forwardscattering_ratio = None
    #     self._data_period = self.parent_sizedist._data_period
    #
    #
    #
    # @property
    # def hemispheric_backscattering(self):
    #     if not self.__hemispheric_backscattering:
    #         out = hemispheric_backscattering(self.angular_scatt_func.data)
    #         out = _timeseries.TimeSeries(out)
    #         out._data_period = self.angular_scatt_func._data_period
    #         self.__hemispheric_backscattering = out
    #     return self.__hemispheric_backscattering
    #
    # @hemispheric_backscattering.setter
    # def hemispheric_backscattering(self,value):
    #     self.__hemispheric_backscattering = value
    #
    # @property
    # def hemispheric_forwardscattering(self):
    #     if not self.__hemispheric_forwardscattering:
    #         out = hemispheric_forwardscattering(self.angular_scatt_func.data)
    #         out = _timeseries.TimeSeries(out)
    #         out._data_period = self.angular_scatt_func._data_period
    #         self.__hemispheric_forwardscattering = out
    #     return self.__hemispheric_forwardscattering
    #
    #
    # @hemispheric_forwardscattering.setter
    # def hemispheric_forwardscattering(self, value):
    #     self.__hemispheric_forwardscattering = value
    #
    # @property
    # def hemispheric_backscattering_ratio(self):
    #     """ratio between backscattering and overall scattering"""
    #     if not self.__hemispheric_backscattering_ratio:
    #         # self.__hemispheric_backscattering_ratio = self.hemispheric_backscattering / self.extinction_coeff
    #         self.__hemispheric_backscattering_ratio = self.hemispheric_backscattering / self.scattering_coeff
    #     return self.__hemispheric_backscattering_ratio
    #
    # @property
    # def hemispheric_forwardscattering_ratio(self):
    #     """ratio between forwardscattering and over scattering"""
    #     if not self.__hemispheric_forwardscattering_ratio:
    #         self.__hemispheric_forwardscattering_ratio = self.hemispheric_forwardscattering / self.scattering_coeff
    #     return self.__hemispheric_forwardscattering_ratio


class OpticalProperties_VP(OpticalProperties):

    @property
    def hemispheric_forwardscattering(self):
        return _vertical_profile.VerticalProfile(super().hemispheric_forwardscattering)

    @property
    def hemispheric_backscattering(self):
        return _vertical_profile.VerticalProfile(super().hemispheric_backscattering)

    @property
    def hemispheric_backscattering_ratio(self):
        return _vertical_profile.VerticalProfile(super().hemispheric_backscattering_ratio)

    @property
    def hemispheric_forwardscattering_ratio(self):
        return _vertical_profile.VerticalProfile(super().hemispheric_forwardscattering_ratio)

    @property
    def absorption_coeff(self):
        return _vertical_profile.VerticalProfile(super().absorption_coeff)

    @property
    def extinction_coeff(self):
        return _vertical_profile.VerticalProfile(super().extinction_coeff)

    @property
    def scattering_coeff(self):
        return _vertical_profile.VerticalProfile(super().scattering_coeff)

    @property
    def _optical_porperties(self):
        if not self._optical_porperties_pv:
            super()._optical_porperties
            layerthickness = self._parent_sizedist.layerbounderies[:, 1] - self._parent_sizedist.layerbounderies[:, 0]
            aod_per_bin_per_layer = self._parent_sizedist.optical_properties.extinction_coeff_per_bin.multiply(layerthickness, axis=0)
            aod_per_layer = _pd.DataFrame(aod_per_bin_per_layer.sum(axis=1), columns=['aod_per_layer'])
            self._aod = aod_per_layer.values.sum()
            aod_cumulative = aod_per_layer.iloc[::-1].cumsum()
            aod_cumulative.rename(columns={'aod_per_layer': 'aod'}, inplace=True)
            self._aod_cumulative = aod_cumulative
        return self._optical_porperties_pv

    @property
    def aod(self):
        self._optical_porperties
        return self._aod

    @property
    def aod_cumulative(self):
        self._optical_porperties
        return _vertical_profile.VerticalProfile(self._aod_cumulative)

class DEPRECATED_OpticalProperties_VP(OpticalProperties):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.extinction_coeff_per_bin = _vertical_profile.VerticalProfile_2D(self.extinction_coeff_per_bin)
        self.aerosol_optical_depth_cumulative_VP = _vertical_profile.VerticalProfile(self._data_dict['AOD_cum'])
        self.asymmetry_param_VP = _vertical_profile.VerticalProfile(self._data_dict['asymmetry_param'])
        self.aerosol_optical_depth_cumulative = self._data_dict['AOD']

class ExtinctionCoeffVerticlProfile(_vertical_profile.VerticalProfile):
    def __init__(self, ext, parent, wavelength, index_of_refraction):
        super(ExtinctionCoeffVerticlProfile, self).__init__(ext)
        self.parent = parent
        self.wavelength = wavelength
        self.index_of_refraction = index_of_refraction

    def plot(self, *args, **kwargs):
        a = super(ExtinctionCoeffVerticlProfile, self).plot(*args, **kwargs)
        a.set_xlabel('Extinction coefficient (m$^{-1}$)')
        return a


def _perform_Miecalculations(diam, wavelength, n, noOfAngles=100.):
    """
    Performs Mie calculations

    Parameters
    ----------
    diam:       NumPy array of floats
                Array of diameters over which to perform Mie calculations; units are um
    wavelength: float
                Wavelength of light in um for which to perform calculations
    n:          complex
                Ensemble complex index of refraction

    Returns
        panda DataTable with the diameters as the index and the mie_scattering results in the different collumns
        total_extinction_coefficient: this takes the sum of all particles crossections of the particular diameter in a qubic
                                      meter. This is in principle the AOD of an L

    """


    diam = _np.asarray(diam)

    extinction_efficiency = _np.zeros(diam.shape)
    scattering_efficiency = _np.zeros(diam.shape)
    absorption_efficiency = _np.zeros(diam.shape)

    extinction_crossection = _np.zeros(diam.shape)
    scattering_crossection = _np.zeros(diam.shape)
    absorption_crossection = _np.zeros(diam.shape)

    # phase_function_natural = pd.DataFrame()
    angular_scattering_natural = _pd.DataFrame()
    # extinction_coefficient = np.zeros(diam.shape)
    # scattering_coefficient = np.zeros(diam.shape)
    # absorption_coefficient = np.zeros(diam.shape)



    # Function for calculating the size parameter for wavelength l and radius r
    sp = lambda r, l: 2. * _np.pi * r / l
    for e, d in enumerate(diam):
        radius = d / 2.

        # print('sp(radius, wavelength)', sp(radius, wavelength))
        # print('n', n)
        # print('d', d)

        mie = _bhmie.bhmie_hagen(sp(radius, wavelength), n, noOfAngles, diameter=d)
        values = mie.return_Values_as_dict()
        extinction_efficiency[e] = values['extinction_efficiency']

        # print("values['extinction_crosssection']",values['extinction_crosssection'])


        scattering_efficiency[e] = values['scattering_efficiency']
        absorption_efficiency[e] = values['extinction_efficiency'] - values['scattering_efficiency']

        extinction_crossection[e] = values['extinction_crosssection']
        scattering_crossection[e] = values['scattering_crosssection']
        absorption_crossection[e] = values['extinction_crosssection'] - values['scattering_crosssection']

        # phase_function_natural[d] = values['phaseFct_natural']['Phase_function_natural'].values
        angular_scattering_natural[float(d)] = mie.get_angular_scatt_func().natural.values

        # print('\n')

    # phase_function_natural.index = values['phaseFct_natural'].index
    angular_scattering_natural.index = mie.get_angular_scatt_func().index

    out = _pd.DataFrame(index=diam)
    out['extinction_efficiency'] = _pd.Series(extinction_efficiency, index=diam)
    out['scattering_efficiency'] = _pd.Series(scattering_efficiency, index=diam)
    out['absorption_efficiency'] = _pd.Series(absorption_efficiency, index=diam)

    out['extinction_crossection'] = _pd.Series(extinction_crossection, index=diam)
    out['scattering_crossection'] = _pd.Series(scattering_crossection, index=diam)
    out['absorption_crossection'] = _pd.Series(absorption_crossection, index=diam)
    return out, angular_scattering_natural


def _get_coefficients(crossection, cn):
    """
    Calculates the extinction, scattering or absorbtion coefficient

    Parameters
    ----------
    crosssection:   float
                    Units are um^2
    cn:             float
                    Particle concentration in cc^-1

    Returns
    --------
    coefficient in m^-1.  This is the differential AOD.
    """
    crossection = crossection.copy()
    cn = cn.copy()
    crossection *= 1e-12  # conversion from um^2 to m^2
    cn *= 1e6  # conversion from cm^-3 to m^-3
    coefficient = cn * crossection

    # print('cn',cn)
    # print('crossection', crossection)
    # print('coeff',coefficient)
    # print('\n')

    return coefficient

def vertical_profile2accumulative_AOD(timeseries):
    data = timeseries.data.copy()
    data.dropna(inplace = True)
    accu_aod = _np.zeros(data.shape)
    for col,k in enumerate(data.keys()):
        series = data[k]
        # series.dropna(inplace = True)
        y = series.values*1e-6
        x = series.index.values

        x = x[::-1]
        y = y[::-1]


        st = 0
        # for e,i in enumerate(x):
        for e in range(x.shape[0]):
            end = e+1
            accu_aod[e][col] = -_integrate.simps(y[st:end], x[st:end])

    accu_aod = _pd.DataFrame(accu_aod, index = x, columns=data.keys())
    accu_aod = _vertical_profile.VerticalProfile(accu_aod)

    accu_aod._x_label = 'AOD$_{abs}$'
    return accu_aod

# from scipy.optimize import curve_fit
from atmPy.aerosols.size_distribution import sizedistribution as sd

class Inversion2SizeDistribution_scenario(object):
    """this handles how ..."""
    def __init__(self, parent, arguments):
        self.parent = parent
        # if arguments == False:
        #     self.size_distribution_parameters = arguments
        
        if isinstance(arguments, _pd.DataFrame):
            assert(False), "is this still used? If you need this, commment out"
            if self.parent.verbose:
                print('arguments are type pd.DatrFrame')
            self.size_distribution_parameters = arguments
            pos = arguments.diameter.values
            numbers = arguments.number.values
            arguments = [val for pair in zip(pos, numbers) for val in pair]
        else:
            if self.parent.verbose:
                print(f'arguments are type {type(arguments)}')
            # if self.parent.verbose:
            #     print('initialize scenario: else')
            # width = self.parent.start_conditions.size_distribution_parameters.width.values # this cased a superloop
            width = self.parent.width_of_aerosol_mode
            self.size_distribution_parameters = _pd.DataFrame({'diameter': arguments[::2], 'number': arguments[1::2], 'width': width})
        
        self._args = arguments
        self.update()


    def update(self):
        # properties
        self._dist = None
        self._extcoeff = None
        self._size_dist_aods = None
        self._size_distribution_moments = None
        self._dist_list = None

    @property
    def args(self):
        return self._args

    @args.setter
    def args(self, value):
        self.update()
        self._args = value

    @property
    def size_distribution(self):
        if isinstance(self._dist, type(None)):
            self._dist = self.args2dist()#self.args)
        return self._dist

    @property
    def size_distribution_list(self):
        if isinstance(self._dist_list, type(None)):
            self.args2dist()
        return self._dist_list

    @property
    def size_distribution_aods(self):
        if isinstance(self._size_dist_aods, type(None)):         
            if isinstance(self.parent.mie_info, type(None)):
                self.AOD_of_simulated_dist()
                
            modes = ('fine', 'coarse')
            aods = _pd.DataFrame(index = modes, columns=self.parent.wavelengths, dtype = _np.float32)
            for wl in self.parent.wavelengths:
                for e,dist in enumerate(self.size_distribution_list):
                    # dist.optical_properties.parameters.refractive_index = self.parent.aerosol_refractive_index
                    # dist.optical_properties.parameters.wavelength = wl
                    dist.optical_properties.parameters.mie_result = self.parent._mie_info[wl]
            
                    aods.loc[modes[e], wl] = dist.optical_properties.extinction_coeff.iloc[0,0]
            self._size_dist_aods = aods
        return self._size_dist_aods

    @property
    def size_distribution_moments(self):
        if isinstance(self._size_distribution_moments, type(None)):
            modes = ('fine', 'coarse')
            moments = ['volume', 'surface']
            moment = _pd.DataFrame(index = modes, columns=moments)
            for e,dist in enumerate(self.size_distribution_list):
                moment.loc[modes[e], 'volume'] = dist.particle_volume_concentration.iloc[0,0]
                moment.loc[modes[e], 'surface'] = dist.particle_surface_concentration.iloc[0,0]
            self._size_distribution_moments = moment
        return self._size_distribution_moments

    def args2dist_deprecated(self, args):
        self.tp_3 = args.copy()
        # print(args)
        assert(args.shape == (4,)), f'shape is {args.shape}'
        positions = args[::2]
        amplitudes = args[1::2]
        dist_list = []
        #     print(args)
        for pos, amp in zip(positions, amplitudes):
            dist = sd.simulate_sizedistribution(
                diameter=self.parent.diameter_range,
                numberOfDiameters=self.parent.number_of_diameters,
                centerOfAerosolMode=pos,
                # widthOfAerosolMode=self.parent.width_of_aerosol_mode,
                widthOfAerosolMode=self.parent.width_of_aerosol_mode,
                numberOfParticsInMode=amp,
            )
            dist_list.append(dist)

        # dist_new = sd.simulate_sizedistribution(new=True)
        # sum em up

        dist = dist_list[0]
        for sdt in dist_list[1:]:
            dist += sdt
        return dist
    
    def args2dist(self):#, args):
        dist_list = []
        
        if self.size_distribution_parameters.iloc[0,0] == _np.nan:
            # create a dummy size dist and set all values to nan
            dist = sd.simulate_sizedistribution(
                                diameter=self.parent.diameter_range,
                                numberOfDiameters=self.parent.number_of_diameters,
                                centerOfAerosolMode=100,
                                # widthOfAerosolMode=self.parent.width_of_aerosol_mode,
                                widthOfAerosolMode=0.2,
                                numberOfParticsInMode=100,
                                )
            dist.data[:] = _np.nan
        
        else:
            for idx, row in self.size_distribution_parameters.iterrows():
                dist = sd.simulate_sizedistribution(
                    diameter=self.parent.diameter_range,
                    numberOfDiameters=self.parent.number_of_diameters,
                    centerOfAerosolMode=row.diameter,
                    # widthOfAerosolMode=self.parent.width_of_aerosol_mode,
                    widthOfAerosolMode=row.width,
                    numberOfParticsInMode=row.number,
                )
                dist_list.append(dist)
    
            # dist_new = sd.simulate_sizedistribution(new=True)
            # sum em up
            self._dist_list = dist_list
            dist = dist_list[0]
            for sdt in dist_list[1:]:
                dist += sdt
        return dist

    @property
    def extinction_coeff(self):
        if isinstance(self._extcoeff, type(None)):
            # self.parent.mie_info
            self._extcoeff = self.AOD_of_simulated_dist()#self.parent.wavelengths, *self.args, )
        return self._extcoeff

    def AOD_of_simulated_dist(self):#, wavelengths, *args, ):

        # args = _np.array(args)
        # dist = self.args2dist()#args)
        dist = self.size_distribution#args)
        channels = self.parent.wavelengths
        first_run_results = {}
        mie_info = self.parent._mie_info
        
        if isinstance(mie_info, type(None)):
            if self.parent.verbose:
                print('doing full mie calculation')

        for wl in channels:
            res = {}
            dt = dist.copy()
            self.parent.tp_dt1 = dt.copy()
            # print(type(self.parent))
            if isinstance(mie_info, type(None)):
                dt.optical_properties.parameters.refractive_index = self.parent.aerosol_refractive_index
                dt.optical_properties.parameters.wavelength = wl
            else:
                dt.optical_properties.parameters.mie_result = mie_info[wl]

            res['extcoeff'] = dt.optical_properties.extinction_coeff
            res['mie_result'] = dt.optical_properties.mie_result.copy()
            # return res
            self.parent.tp_dt2 = dt.copy()
            first_run_results[wl] = res

        if isinstance(mie_info, type(None)): 
            self.parent._mie_info = {i: v['mie_result'] for i, v in first_run_results.items()}

        ext_coeff = _np.array([first_run_results[i]['extcoeff'].iloc[0, 0] for i in first_run_results])
        return ext_coeff



class Inversion2SizeDistribution(object):
    def __init__(self, sfr_AOD_test_dp, # this is a row of a pandas frame
                 diameter_range = [80, 20000], 
                 number_of_diameters = 200,
                 aerosol_refractive_index = 1.5,
                 width_of_aerosol_mode = (0.2, 0.25),
                 # start_conditions = [400, 2000., 800, 320],
                 start_args = [120, 2000.0, 1000, 1],
                 pre_fit_amps = True,
                 # bounds = [(-_np.inf, -_np.inf), (-_np.inf, -_np.inf), (-_np.inf, -_np.inf), (-_np.inf, -_np.inf)],
                 bounds=[(10, 0, 400, 0), (400, _np.inf, 6000, _np.inf)],
                 mie_info = None,
                 verbose = False):
        """
        

        Parameters
        ----------
        sfr_AOD_test_dp : TYPE
            DESCRIPTION.
        # this is a row of a pandas frame                 diameter_range : TYPE, optional
            DESCRIPTION. The default is [80, 30000].
        number_of_diameters : TYPE, optional
            DESCRIPTION. The default is 100.
        aerosol_refractive_index : TYPE, optional
            DESCRIPTION. The default is 1.5.
        width_of_aerosol_mode : float or array-like, optional
            If single number it is applied to all aerosol modes. To apply individual width use array-like type. The default is 0.15.
        # start_conditions : TYPE, optional
            DESCRIPTION. The default is [400, 2000., 800, 320].
        start_args : TYPE, optional
            DESCRIPTION. The default is [150, 2000.0, 1500, 1].
        # bounds : TYPE, optional
            DESCRIPTION. The default is [(-_np.inf, -_np.inf), (-_np.inf, -_np.inf), (-_np.inf, -_np.inf), (-_np.inf, -_np.inf)].
        bounds : TYPE, optional
            DESCRIPTION. The default is [(10, 0, 800, 0), (800, _np.inf, 6000, _np.inf)].
        mie_info : TYPE, optional
            DESCRIPTION. The default is None.
        verbose : TYPE, optional
            DESCRIPTION. The default is False.

        Returns
        -------
        None.

        """
        self.diameter_range =diameter_range
        self.number_of_diameters = number_of_diameters
        self.aerosol_refractive_index = aerosol_refractive_index
        self.width_of_aerosol_mode = width_of_aerosol_mode
        self.sfr_AOD_test_dp = sfr_AOD_test_dp.copy().astype(_np.float64)
        self.wavelengths = self.sfr_AOD_test_dp.index
        self.verbose = verbose
        self.pre_fit_amps = pre_fit_amps
        self.fit_bounds = _np.array(bounds)
        self.fit_cutoff = 1e-10
        self.fit_cutoff_xtol = None #self.fit_cutoff
        self.fit_cutoff_ftol = None #self.fit_cutoff
        self.fit_cutoff_gtol = self.fit_cutoff
        self.fit_scale_vars_amp = [1,1]
        self.fit_scale_vars = [1,1,1,1]
        self.fit_max_nfev = 20
        self.fit_tr_solver = None
        self.fit_jac = '2-point'
        self.fit_method = 'dogbox' #trf
        self.fit_diff_step = None
        self.rerun_if_fail = True #if res is not good enough trf will be run untill convergence

        # fit_amp settings
        self.fit_amp_problem_type = 'lin' # lin, nolin (linear or not linear), lin is 5 times faster

        #properties
        self._start_conditions = None
        # self.start_conditions = start_conditions
        self.start_args = start_args
        self._fit_result = None
        self._mie_info = mie_info
        self.scale = 'log' # if optimization is done in log scale, should be preferable! Parameters and extinction will be handled as log scale!!
        


    @property
    def mie_info(self):
        """
        Performs Mie calculations for all diameter steps, so it is not calcualted at each iteration.

        Returns
        -------
        TYPE
            DESCRIPTION.

        """
        # if isinstance(self._mie_info, type(None)):
        #     if self.verbose:
        #         print('getting mieinfo')
        #     out = self.start_conditions.AOD_of_simulated_dist(self.sfr_AOD_test_dp.index, *self.start_conditions.args, is_initial_run=True)
        #     # self._mie_info = {i: v['mie_result'] for i, v in out['first_run_results'].items()}
        #     # self.start_conditions.args = out['args']
        #     if self.verbose:
        #         print('done mieinfo')
        return self._mie_info

    @property
    def start_conditions(self):
        if isinstance(self._start_conditions, type(None)):
            # print('hasdfasdfasd')
            # args = [400, 2000., 800, 320] # if start conditions are not set this will be used ad a default
            # self._start_conditions = Inversion2SizeDistribution_scenario(self, self.start_args)
            self.start_conditions = self.start_args #note, this will execute the setter!!!
        return self._start_conditions

    @start_conditions.setter
    def start_conditions(self, value):
        """
        set the start conditions

        Parameters
        ----------
        value : list [pos1, particle_nu_1, pos2, part.....]
            the default is [400, 2000., 800, 320].

        Returns
        -------
        None.

        """
        if self.verbose:
            print('setting start_conditions')
        if isinstance(value, type(None)):
            if self.verbose:
                print('start condition value is None, use default')
            arguments = self.start_args.copy()#[400, 2000., 800, 320] # if start conditions are not set this will be used ad a default
        else:
            arguments = value # if start conditions are not set this will be used ad a default
        pre_start_conditions = Inversion2SizeDistribution_scenario(self, arguments)
        
        # We need to scale the particle numbers such that we are already close to observed values. Otherwise the optimization will run away.
        assert(self.fit_amp_problem_type in ['lin', 'nolin'])
        if self.fit_amp_problem_type == 'nolin': #not needed if i do the linearization of the fit_amp.
            factor = self.sfr_AOD_test_dp.values[1] / pre_start_conditions.extinction_coeff[1]# adjust values on 500 channel only
            arguments = _np.array(arguments, dtype = float)
            arguments[1::2] *= float(factor)
        # else:
        #     self.size_distribution_parameters = _pd.DataFrame({'diameter': arguments[::2], 'number': arguments[1::2], 'width': width})
        #     self.size_distribution_parameters = _pd.DataFrame({'diameter': arguments[::2], 'number': arguments[1::2], 'width': width})
        
        self._start_conditions = Inversion2SizeDistribution_scenario(self, arguments)

    @property
    def fit_result(self):
        if isinstance(self._fit_result, type(None)):
            start = _pd.Timestamp.now()
            fitrs = self.fit()
            if fitrs == False:
                args = [_np.nan]*4
                self._fit_result = Inversion2SizeDistribution_scenario(self, args)
                self._fit_result.sigma = _np.nan
            else:
                fitrsx = fitrs.x
                if self.scale == 'log':
                    # fitrsx = 10**fitrsx
                    fitrsx[::2] = _np.power(10, fitrsx[::2])
                fit_res = Inversion2SizeDistribution_scenario(self, fitrsx)
                fit_res.full_result = fitrs
                fit_res.cpu_usage = fitrs.cpu_usage
                fit_res.sigma = _np.sqrt(((self.sfr_AOD_test_dp.values-_np.array(fit_res.extinction_coeff))**2).sum())
                self._fit_result = fit_res
            
            self._fit_result.optimization_time = _pd.Timestamp.now() - start
        return self._fit_result

    def cost_fun(self, *params):
        measured = self.sfr_AOD_test_dp.values
        args = params[0].copy()
        self.tp_params = args
        # assert(False)
        if self.scale == 'log':
            # params = _np.array(params)
            
            # self.tp_params = params
            # assert(False)
            # params = 10**params
            # params = _np.power(10,params)
            args[::2] = _np.power(10, args[::2])
            self.tp_params = args
            # assert(False)
        self.tp_params = params
        # assert(False)

        scene = Inversion2SizeDistribution_scenario(self, args)
        ext = scene.extinction_coeff
        
        # if self.scale == 'log':
        if 1:
            measured = _np.log10(measured)
            ext = _np.log10(ext)
        
        cost = measured - ext
        self.tp_cost.append((cost**2).sum())
        self.tp_args.append(scene.args)
        #### create a scalar cost value
        # for some reason those scalar vales are not creating a good result. I suggest to stay with the 1d array
        # cost = float(_np.sum(cost)) # does not work well
        # cost = float(_np.mean(cost))
        # cost = float(_np.std(cost))
        self.no_of_iterations += 1
        return cost
    
    def fit(self):
        if _np.isnan(self.sfr_AOD_test_dp.values).sum() > 0:
            # out = _np.zeros(self.start_conditions.args.shape[0])
            # out = _np.zeros(self.start_conditions.args.shape[0])
            # out[:] = _np.nan
            return False
        
        #### computations cost start
        self.no_of_iterations = 0
        p = psutil.Process(os.getpid())
        start_time = time.time()
        p.cpu_percent()
        ###
        
        self.tp_cost = []
        self.tp_args = []
        
        if self.pre_fit_amps:
            arguments = self.fit_amps().args.copy()
        else:
            assert(False), 'This is not recommended and hardcoded to True'
            arguments = self.start_conditions.args
            
        
        bounds = self.fit_bounds.copy()
        if self.scale == 'log':
            # arguments = _np.log10(arguments)
            # bounds = _np.log10(bounds)
            
            arguments[::2] = _np.log10(arguments[::2])
            def replwl(b):
                b[::2] = _np.log10(b[::2])
                return b
            bounds = _np.array([replwl(b) for b in bounds])
        
        # if isinstance(self.fit_scale_vars, type(None)):
        #     self.fit_scale_vars = arguments
            
        # if isinstance(self.fit_diff_step , type(None)):
        #     self.fit_diff_step = _np.array([5,1,2,1]) * 0.00005
            
        self.tp_arguments = arguments
        self.tp_bounds = bounds
        try:
            fit_res = scipy.optimize.least_squares(self.cost_fun,
                                        # self.start_conditions.args,
                                        arguments,
                                    #                     sigma=None,
                                        #                     absolute_sigma=False,
                                        #                     check_finite=True,
                                        bounds=bounds,
                                        method=self.fit_method,
                                        jac = self.fit_jac,
                                        xtol=self.fit_cutoff_xtol,
                                        ftol=self.fit_cutoff_ftol,
                                        gtol=self.fit_cutoff_gtol,
                                        verbose=0,
                                        x_scale = self.fit_scale_vars,
                                        max_nfev = self.fit_max_nfev,
                                        tr_solver = self.fit_tr_solver,
                                        diff_step= self.fit_diff_step,  
                                        # this merely changes the initial steps and its important only if the steps are to big and jump over minima, not the case here
                                        # Older comment: Although it sounds like its something different... it works!!!!!!
                                        # I don't think this is what I thought it is!! This has something to do with the computers resolution
                                        # Removed it none the less. I think the initial adjustment of amplitudes makes this unnecessary and I see better results witout it.
                                        #                     col_deriv=True
                                        )
        except ValueError as e:
            print(str(e))
            fit_res = type('adhoc_fitres', (), {'success': False,
                                                'x': arguments,
                                                'cost': _np.nan})
            
        if not fit_res.success:
            if self.fit_method == 'dogbox' and self.rerun_if_fail:
                try:
                    fit_res = scipy.optimize.least_squares(self.cost_fun,
                                            fit_res.x,
                                            bounds=bounds,                                            
                                            jac = self.fit_jac,
                                            method='trf',
                                            xtol=self.fit_cutoff,
                                            ftol=self.fit_cutoff,
                                            gtol=self.fit_cutoff,
                                            x_scale = self.fit_scale_vars,
                                            # max_nfev = self.max_nfev,
                                            )
                except ValueError as e:
                    print(str(e))
                    fit_res = type('adhoc_fitres', (), {'success': False,
                                                        'x': arguments,
                                                        'cost': _np.nan})
        self.tp_fit_res = fit_res
        #### computations cost finallize
        cpu_percent = p.cpu_percent()
        end_time = time.time()
        execution_time = end_time - start_time
        accumulated_cpu_usage = cpu_percent * execution_time
        fit_res.cpu_usage = accumulated_cpu_usage
        return fit_res

    def cost_fun_amps(self, amps):
        """
        Deprecated: this should not be necessary anymore. It is needed to do a non-linear optimization, 
        which was replaced by a linear approach. This is kept for the eventuallity that the linear approach shows weaknesses.
        This fits only the amplitudes but not the positions ... just to get started!.

        Parameters
        ----------
        *params : TYPE
            DESCRIPTION.

        Returns
        -------
        cost : TYPE
            DESCRIPTION.

        """
        params = self.start_args.copy()
        # print(amps)
        # print(type(amps))
        params[1::2] = amps
        measured = self.sfr_AOD_test_dp.values
        scene = Inversion2SizeDistribution_scenario(self, params).extinction_coeff
        # measured = _np.log10(measured)
        # scene = _np.log10(scene)
        cost = measured - scene
        # self.tp_cost.append(_np.mean(cost))
        #### create a scalar cost value
        # for some reason those scalar vales are not creating a good result. I suggest to stay with the 1d array
        # cost = float(_np.sum(cost)) # does not work well
        # cost = float(_np.mean(cost))
        # cost = float(_np.std(cost))
        self.tp_fit_amp_iterations += 1
        return cost
    
    def fit_amps(self):
                    
        arguments = self.start_conditions.args.copy()
        
        if _np.isnan(self.sfr_AOD_test_dp.values).sum() > 0:
            # out = _np.zeros(self.start_conditions.args.shape[0])
            # out = _np.zeros(self.start_conditions.args.shape[0])
            # out[:] = _np.nan
            return False
        
        # self.tp_cost = []
        # ds = 
        # mode_scale = 10
        # ds = 0.51
        # scale_vars = True
        
        # if self.scale_vars:
        #     amp_f = arguments[1]
        #     amp_c = arguments[3]
        #     costf = lambda amp_f, amp_c: (self.cost_fun_amps([amp_f, amp_c])**2).sum()**(1/2) # this is not exactly the same cost as what least_sqrt does
            
        #     cost_all = costf(amp_f, amp_c)
        #     cost_f = abs(costf(amp_f+10, amp_c) - cost_all)
        #     cost_c = abs(costf(amp_f, amp_c+10) - cost_all)
        #     scale = [1, cost_f / cost_c]
        # else:
        #     scale = [1,1]
        
        self.tp_fit_amp_iterations = 0
        scale = self.fit_scale_vars_amp
        # print(scale)
        params = self.start_args.copy()
        if self.fit_amp_problem_type == 'nolin':
            fit_res = scipy.optimize.least_squares(self.cost_fun_amps,
                                            # self.start_conditions.args,
                                            arguments[1::2],
                                        #                     sigma=None,
                                            #                     absolute_sigma=False,
                                            #                     check_finite=True,
                                            # bounds=self.bounds,
                                            method='trf',
                                                                # method='lm',
                        #                                         jac=None,
                                            #                     maxfev=1,
                                            xtol=self.fit_cutoff,
                                            ftol=self.fit_cutoff,
                                            gtol=self.fit_cutoff,
                                            #                     callback=lambda x: print(x)
                                            #                     dict(kwx_scale = 'jac'),
                                            verbose=0,
                                            x_scale = scale,
                                            # x_scale='jac',
                                            # x_scale=[1,0.2,1,0.2],
                                            # x_scale=_np.array([1,1.5,1, mode_scale]),#*1e-3,
                                            # diff_step=[ds ** 4, ds, ds ** 4, ds] # Although it sounds like its something different... it works!!!!!!I don't think this is what I thought it is!! This has something to do with the computers resolution
                                            #                     col_deriv=True
                                            )
            params[1::2] = fit_res.x
        elif self.fit_amp_problem_type == 'lin':
            df = self.start_conditions.size_distribution_aods.transpose()
            df['target'] = self.sfr_AOD_test_dp
            
            A = df.loc[:,['fine', 'coarse']].values.astype(float)
            b = df.target.values
            
            fit_res = scipy.optimize.lsq_linear(A, b, bounds = (0, _np.inf))
            params[1::2] *= fit_res.x
                     
        else:
            assert(self.fit_amp_problem_type in ['lin', 'nolin'])
            
        
        self.tp_fitres_amp = fit_res
        fit_result = Inversion2SizeDistribution_scenario(self, params)
        fit_result.full_result = fit_res
        # fit_result.sigma = _np.sqrt(((self.sfr_AOD_test_dp.values-_np.array(self.fit_result.extinction_coeff))**2).sum())
        return fit_result


