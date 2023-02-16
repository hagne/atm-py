# import pkg_resources
import xarray as _xr
# import importlib.resources
# import pkgutil
# from io import StringIO
# import os
import pathlib as pl
import numpy as np
import xarray as xr
import scipy as sp

class EquivalentDiameters(object):
    """ Class that handles converstion of various equivalent particle diameters."""
    def __init__(self,
                 diameter = None,
                 diameter_equivalent = None, 
                 density_measured = None, 
                 density_calibrated = None, 
                 dynamic_shape_factor = None, 
                 shape_parameter = None, 
                 axis_ratio = None, 
                 ):
        self._reset()
        if diameter_equivalent == 'volume':
            self._dvol = diameter
        elif diameter_equivalent == 'aerodynamic':
            self._daero = diameter
        elif diameter_equivalent == 'mobility':
            self._dmobil = diameter
        elif isinstance(diameter_equivalent, type(None)):
            if not isinstance(diameter, type(None)):
                self._dvol = diameter
        else:
            raise ValueError(f'{diameter_equivalent} is not a valid option, chose from: volume, aerodynamic, mobility.')
          
            
        # self.diameter_equivalent =diameter_equivalent
        self.density_measured = density_measured
        self.density_calibrated = density_calibrated
        self.dynamic_shape_factor = dynamic_shape_factor 
        self.shape_parameter = shape_parameter
        self.axis_ratio = axis_ratio
        
        assert(sum([not isinstance(b, type(None)) for b in [dynamic_shape_factor, shape_parameter, axis_ratio]]) <=1), 'Only one of the particle shape discribing argurments (dynamic_shape_factor, shape_parameter, axis_ratio) can be set.'
    
    def _reset(self):
        self._dvol = None
        self._daero = None
        self._dmobil = None
    
    @property
    def volume_diameter(self):
        if isinstance(self._dvol, type(None)):
            if not isinstance(self._daero, type(None)):
                self._dvol = self._daero * self._d_aps2d_vol()
            elif not isinstance(self._dmobil, type(None)):
                self._dvol = self._dmobil * self._dmobil2dvol()
        return self._dvol
    
    @volume_diameter.setter
    def volume_diameter(self, value):
        self._reset()
        self._dvol = value
        return
    
    @property
    def aerodynamic_diameter(self):
        if isinstance(self._daero, type(None)):
            if not isinstance(self._dvol, type(None)):
                self._daero = self._dvol / self._d_aps2d_vol()
            elif not isinstance(self._dmobil, type(None)):
                assert(False), 'not implemented yet!!'
        return self._daero
    
    @aerodynamic_diameter.setter
    def aerodynamic_diameter(self, value):
        self._reset()
        self._daero = value
        return
    
    @property
    def mobility_diameter(self):
        if isinstance(self._dmobil, type(None)):
            if not isinstance(self._daero, type(None)):
                assert(False), 'not implemented yet!!'
            elif not isinstance(self._dvol, type(None)):
                assert(False), 'not implemented yet!!'
        return self._dmobil
    
    @mobility_diameter.setter
    def mobility_diameter(self, value):
        self._reset()
        self._dmobil = value
        return
    
    def _dmobil2dvol(self):
        xsi = self.dynamic_shape_factor
        assert(not isinstance(xsi, type(None))), 'dynamic_shape_factor needs to be set.'
        return 1/xsi
    
    def _d_aps2d_vol(self):#, roh_p = 2, roh_0 = 2, xsi = 1):
        """
        Provides a correction factor to convert the aerodynamic diameter to the volume equivalent diameter
        
        Parameters
        ----------
        roh_p: actual particles density
        roh_0: calibration density (either the acuatl calibration material, e.g. PSL, but more often then not an additional adjustem to a more realistic density, e.g. 2 in case of SGP
        xsi: dynamic shape factor
        
        Returns
        -------
        correction factor for volume equivalent diameter
        """
        roh_0 = self.density_calibrated
        roh_p = self.density_measured
        xsi = self.dynamic_shape_factor
        assert(not isinstance(roh_0, type(None))), 'density_calibrated needs to be set.'
        assert(not isinstance(roh_p, type(None))), 'density_measured needs to be set.'
        assert(not isinstance(xsi, type(None))), 'dynamic_shape_factor needs to be set.'
        c = np.sqrt((roh_0 * xsi) / roh_p)
        return c

def _load_davies1979():
    """
    This dataset is an extraction of data from 
    
    Davies, C. N. (1979). Particle-fluid interaction. Journal of Aerosol 
    Science, 10(5), 477–513. https://doi.org/10.1016/0021-8502(79)90006-5
    
    For details on how the dataset extracted/generated check the notebook
    fundgrube/aerosols/sphericity/aspectratio_vs_dynamicshapefactor.ipynb

    Returns
    -------
    out : TYPE
        DESCRIPTION.

    """
    p2f = pl.Path(__file__).parent.joinpath('data/spheroid2dynamic_shape_Davies1979.nc')
    out = _xr.open_dataset(p2f)
    return out

class AspectRatio2DynamicShapeFactor():
    def __init__(self,basedon="davies1979"):
        self.basedon = basedon
        self._dataset = None
        return
    
    @property
    def dataset(self):
        if isinstance(self._dataset, type(None)):
            assert(self.basedon == 'davies1979')
            self._dataset = _load_davies1979()
        return self._dataset
    
    def shape_distribution2effective_dsf(self, limits, distribution_function = 'squared', 
                                         # scale = 5,
                                         ):
        """
        Determines the effective dynamic shape factor from a distribution of
        aspect ratios based on the shape parameter.

        Parameters
        ----------
        limits : TYPE
            DESCRIPTION.
        distribution_function : TYPE, optional
            DESCRIPTION. The default is 'squared'.

        Returns
        -------
        dsf_eff : TYPE
            DESCRIPTION.

        """
        # testing if sampling rate is uniform
        s2d = self.dataset
        nda = s2d.shape_parameter.values
        spacing = nda[1:] - nda[:-1]
        assert(np.std(spacing)/np.mean(spacing) < 1e-3), 'Data is not equally distributed, this will case a bias!! You probably changed the data at some point, because that did not happen with the original data!!!'
        # limits
        # if limits array is only one element long (happens in optimization routine)
        if hasattr(limits, '__iter__'):
            if len(limits) < 2:
                limits = limits[0]
         
        # if limits is single number make limits symetric around 0
        if not hasattr(limits, '__iter__'):
            limits = [-limits,limits]
        
        
        assert(limits[0] > float(self.dataset.shape_parameter.min())), 'limits[0] out of bounds'
        assert(limits[1] < float(self.dataset.shape_parameter.max())), 'limits[1] out of bounds'
        self.tp_limits = limits
        
            
        
        # effective dynamicshapefactor
        if distribution_function == 'squared':
            # nofv = (limits[1] - limits[0])/np.mean(spacing) * scale
            #int(np.ceil(nofv)) 
            # self.tp_nsp = int(np.ceil(nofv))
            
            noofv = 250 # the number of 250 results in an error smaller than 0.5% for dsf between 1.05, and 1.5
            new_shape_params = np.linspace(limits[0], limits[1], noofv)
            new_dsf = self.dataset.dynamic_shape_factor.interp(shape_parameter=new_shape_params, method = 'quadratic')
            self.tp_sp_dist = new_dsf
            dsf_eff = float(new_dsf.sum()/new_dsf.shape[0])
            
        elif distribution_function == 'squared1':
            sp_dist = np.ones(s2d.shape_parameter.shape)
            sp_dist[np.where(s2d.shape_parameter < limits[0])] = 0
            sp_dist[np.where(s2d.shape_parameter > limits[1])] = 0
            
            sp_dist = xr.DataArray(sp_dist, 
                    coords=[s2d.shape_parameter,],
                    )
            
            self.tp_sp_dist = sp_dist
            dsf_eff = float((s2d.dynamic_shape_factor * sp_dist).sum()/ sp_dist.sum())
            
        else:
            assert(False), 'not implemented yet'
            
        return dsf_eff
            
    def effective_dsf2shape_distribution(self, dsf, bounds=(1, 10)):
        res = sp.optimize.least_squares(lambda sp: (self.shape_distribution2effective_dsf(sp) - dsf) * 1e2, [5,], bounds=bounds) 
        return float(res.x)

# Sphereoids
def volume_equivalent_radius_spheroid(equatorial_axis, polar_axis):
    a = equatorial_axis
    c = polar_axis
    return  (a**2 * c)**(1/3)

def projection_area_equivalent_radius_spheroid(axis1, axis2):
    a = axis1
    c = axis2
    return (a * c) ** (1 / 2)


