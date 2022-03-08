from copy import deepcopy

import numpy as _np
import matplotlib.pylab as _plt
from matplotlib.colors import LogNorm
from matplotlib.ticker import MaxNLocator as _MaxNLocator

from atmPy.tools import plt_tools, math_functions, array_tools
from atmPy.tools import pandas_tools as _panda_tools
from atmPy.general import timeseries as _timeseries
from atmPy.general import vertical_profile as _vertical_profile
import pandas as pd
import warnings as _warnings
import datetime
import scipy as _sp
import scipy.optimize as optimization
from scipy import stats
from atmPy.aerosols.physics import hygroscopicity as hygroscopicity
from atmPy.tools import pandas_tools
from atmPy.aerosols.physics import optical_properties
from atmPy.aerosols.size_distribution import moments
from atmPy.gases import physics as _gas_physics
from . import modes
from atmPy.aerosols.size_distribution import diameter_binning as _diameter_binning
try:
    from statsmodels import robust as _robust
except ModuleNotFoundError:
    _warnings.warn('statsmodels not installed. You might encounter some functionality limitations.')
try:
    import xarray as _xr
except ModuleNotFoundError:
    _warnings.warn('xarray not installed. You might encounter some functionality limitations.')
# from matplotlib.ticker import MultipleLocator
from matplotlib import gridspec
# from atmPy import atmosphere


# Todo: rotate the plots of the layerseries (e.g. plot_particle_concentration) to have the altitude as the y-axes

# # TODO: Fix distrTypes so they are consistent with our understanding.
# distTypes = {'log normal': ['dNdlogDp', 'dSdlogDp', 'dVdlogDp'],
#              'natural': ['dNdDp', 'dSdDp', 'dVdDp'],
#              'number': ['dNdlogDp', 'dNdDp'],
#              'surface': ['dSdlogDp', 'dSdDp'],
#              'volume': ['dVdlogDp', 'dVdDp']}

_axes_types = ('AxesSubplot', 'AxesHostAxes')

_colors = _plt.rcParams['axes.prop_cycle'].by_key()['color']

def align2sizedist(sizedist, other):
    if type(other).__name__ in ('int', 'float'):
        pass
    elif type(other).__name__  in ('TimeSeries'):
        if not _np.array_equal(sizedist.data.index, other.data.index):
            other = other.align_to(sizedist)
        other = other.data

    if type(other).__name__ in ('DataFrame', 'ndarray'):
        if other.shape[0] != sizedist.data.shape[0]:
            txt = """\
    Length of new array has to be the same as that of the size distribution. Use
    sizedistribution.align to align the new array to the appropriate length"""
            raise ValueError(txt)

        if not _np.array_equal(sizedist.data.index, other.index):
            txt = """\
    The index of the new DataFrame has to be the same as that of the size distribution. Use
    sizedistribution.align to align the index of the new array."""
            raise ValueError(txt)
    return other

def merge_size_distributions(dist_self, dist_other, fill_value = 0, round_dec = 5):
    """
    Experimental!!! Merges (adds) two sizedistributions that have different length. Currently this is only working if the
    overlapping section of the size distributions is aligned (bins are exactly the same where overlapping)

    Parameters
    ----------
    dist_self
    dist_other
    fill_value: float [0]
        When adding a value to nan the result is nan. Therefore all the nans are replaced with this value.
    round_dec: int [5]
        sometimes the bins are not equal due to rounding issue in the 10 or so digit. Rounding to the 5th digit ususally
        takes care of that without introducing an error

    Returns
    -------
    size distribution

    """
    assert(dist_self.distributionType == dist_other.distributionType)

    dist_self.bins = _np.round(dist_self.bins, round_dec)
    dist_other.bins = _np.round(dist_other.bins, round_dec)
    dist_self.data[_np.isnan(dist_self.data)] = 0
    dist_other.data[_np.isnan(dist_other.data)] = 0
    data =pd.DataFrame.add(dist_self.data, dist_other.data, fill_value=fill_value)
    bins_merged = _np.unique(_np.sort(_np.append(dist_self.bins, dist_other.bins)))
    dist_out = SizeDist(data, bins_merged, dist_self.distributionType)
    return dist_out

def fit_normal_distribution2sizedist(sizedist, log=True, p0=[10, 180, 0.2], show_error = False, curve_fit_kwargs = None):
    """ Fits a single normal distribution to each line in the data frame.

    Parameters
    ----------
    p0: array-like
        fit initiation parameters [amp, pos, width(log-width)]
    curve_fit_kwargs: dict
        Additional kwargs that are  passed to the fit routine

    log: not really working

    Returns
    -------
    pandas DataFrame instance (also added to namespace as data_fit_normal)

    """

    def fit_normal_dist(x, y, log=True, p0=[10, 180, 0.2], test=False, curve_fit_kwargs=None):
        """Fits a normal distribution to a """

        if type(curve_fit_kwargs) == type(None):
            curve_fit_kwargs = {}

        param = p0[:]
        x = x[~ _np.isnan(y)]
        y = y[~ _np.isnan(y)]
        if x.shape == (0,):
            return [_np.nan] * 5

        if log:
            x = _np.log10(x)
            param[1] = _np.log10(param[1])
            if 'bounds' in curve_fit_kwargs.keys():
                curve_fit_kwargs['bounds'][0][1] = _np.log10(curve_fit_kwargs['bounds'][0][1])
                curve_fit_kwargs['bounds'][1][1] = _np.log10(curve_fit_kwargs['bounds'][1][1])

        # try:
        para = optimization.curve_fit(math_functions.gauss, x.astype(float), y.astype(float), p0=param,
                                      **curve_fit_kwargs)
        # except:
        #     para = optimization.curve_fit(math_functions.gauss, x.astype(float), y.astype(float), p0=param,
        #                                   **curve_fit_kwargs)

        amp = para[0][0]
        sigma = para[0][2]
        if log:
            pos = 10 ** para[0][1]
            sigma_high = 10 ** (para[0][1] + para[0][2])
            sigma_low = 10 ** (para[0][1] - para[0][2])
        else:
            pos = para[0][1]
            sigma_high = (para[0][1] + para[0][2])
            sigma_low = (para[0][1] - para[0][2])
        if test:
            return [amp, pos, sigma, sigma_high, sigma_low], para
        else:
            return [amp, pos, sigma, sigma_high, sigma_low]

    if type(curve_fit_kwargs) == type(None):
        curve_fit_kwargs = {}

    sd = sizedist.copy()

    if sd.distributionType != 'dNdlogDp':
        if sd.distributionType == 'calibration':
            pass
        else:
            _warnings.warn(
                "Size distribution is not in 'dNdlogDp'. I temporarily converted the distribution to conduct the fitting. If that is not what you want, change the code!")
            sd = sd.convert2dNdlogDp()

    n_lines = sd.data.shape[0]
    amp = _np.zeros(n_lines)
    pos = _np.zeros(n_lines)
    sigma = _np.zeros(n_lines)
    sigma_high = _np.zeros(n_lines)
    sigma_low = _np.zeros(n_lines)
    for e, lay in enumerate(sd.data.values):
        if show_error:
            fit_res = fit_normal_dist(sd.bincenters, lay, log=log, p0=p0, curve_fit_kwargs=curve_fit_kwargs)
        else:
            try:
                fit_res = fit_normal_dist(sd.bincenters, lay, log=log, p0=p0, curve_fit_kwargs = curve_fit_kwargs)
            except (ValueError, RuntimeError):
                fit_res = [_np.nan, _np.nan, _np.nan, _np.nan, _np.nan]
        amp[e] = fit_res[0]
        pos[e] = fit_res[1]
        sigma[e] = fit_res[2]
        sigma_high[e] = fit_res[3]
        sigma_low[e] = fit_res[4]

    df = pd.DataFrame()
    df['Amp'] = pd.Series(amp)
    df['Pos'] = pd.Series(pos)
    df['Sigma'] = pd.Series(sigma)
    df['Sigma_high'] = pd.Series(sigma_high)
    df['Sigma_low'] = pd.Series(sigma_low)
    df.index = sd.data.index
    # sizedist.data_fit_normal = df
    if type(sizedist).__name__ == 'SizeDist_TS':
        out = _NormalDistributionFitRes_TS(sizedist, df, sampling_period=sizedist._data_period)
    elif type(sizedist).__name__ == 'SizeDist_LS':
        out = _NormalDistributionFitRes_VP(sizedist, df)
    else:
        out = df
    return out

def open_csv(fname, fill_data_gaps_with = None, ignore_data_gap_error = False,
             # fixGaps=False
             ):
    """

    Args:
        fname:
        fill_data_gaps_with: float
            If None gaps are not filled. This should eighter be np.nan (if the instrument failed) or 0 if particle
            concentration was so low, that no particle was detected in that time window
        ignore_data_gap_error:

    Returns:

    """
    headerNo = 50
    rein = open(fname, 'r')
    nol = ['distributionType', 'objectType']
    outDict = {}
    for i in range(headerNo):
        split = rein.readline().split('=')
        variable = split[0].strip()
        if split[0][0] == '#':
            break
        value = split[1].strip()
        if variable in nol:
            outDict[variable] = value
        else:
            outDict[variable] = _np.array(eval(value))
        if i == headerNo - 1:
            raise TypeError('Sure this is a size distribution?')

    rein.close()
    data = pd.read_csv(fname, header=i + 1, index_col=0)
    # data.index = pd.to_datetime(data.index)
    if outDict['objectType'] == 'SizeDist_TS':
        data.index = pd.to_datetime(data.index)
        distRein = SizeDist_TS(data, outDict['bins'], outDict['distributionType'],
                               fill_data_gaps_with=fill_data_gaps_with,
                               ignore_data_gap_error=ignore_data_gap_error
                               # fixGaps=fixGaps
                               )
    elif outDict['objectType'] == 'SizeDist':
        distRein = SizeDist(data, outDict['bins'], outDict['distributionType'],
                            # fixGaps=fixGaps
                            )
    elif outDict['objectType'] == 'SizeDist_LS':
        distRein = SizeDist_LS(data, outDict['bins'], outDict['distributionType'],
                               # fixGaps=fixGaps
                               )
    else:
        raise TypeError('not a valid object type')
    return distRein


def open_netcdf(fname):
    ds = _xr.open_dataset(fname, autoclose=True)
    data = ds.sizedistribution.to_pandas()
    binedges = ds.binedges.data
    moment = ds.sizedistribution.moment
    disttype = ds.attrs['type']# ds._atmPy.to_pandas().loc['type']

    if 'normal_distribution_fits' in ds.variables.keys():
        normfitres = ds.normal_distribution_fits.to_pandas()
    else:
        normfitres = None

    if 'housekeeping' in ds.variables.keys():
        hk = ds.housekeeping.to_pandas()
    else:
        hk = None

    if  disttype == 'SizeDist_TS':
        dist = SizeDist_TS(data, binedges, moment, ignore_data_gap_error=True)
        dist._data_period = float(ds.attrs['data_period'])#  ds._atmPy.loc['data_period'].values)
        if type(hk) != type(None):
            dist.housekeeping = _timeseries.TimeSeries(hk,
                                                       sampling_period=dist._data_period)
        if type(normfitres) != type(None):
            dist.normal_distribution_fits = _NormalDistributionFitRes_TS(dist, normfitres,
                                                                   sampling_period=dist._data_period)

    elif  disttype == 'SizeDist_LS':
        dist = SizeDist_LS(data, binedges, moment, ds.layer_boundaries.to_pandas().values)
        if type(hk) != type(None):
            dist.housekeeping = _vertical_profile.VerticalProfile(hk)
        if type(normfitres) != type(None):
            dist.normal_distribution_fits = _NormalDistributionFitRes_VP(dist, normfitres)

    else:
        raise ValueError('open_netcdf is not defined for type {} yet ... programming requried'.format(disttype))

    ## corrections
    corrections = ds.applied_corrections.split(', ')

    if 'flow_in_LFE4temp_difference' in corrections:
        dist._correct4ambient_LFE_tmp_difference = True
    if 'flow_rate' in corrections:
        dist._correct4fowrate = True

    return dist

def get_label(distType):
    """ Return the appropriate label for a particular distribution type
    """
    if distType == 'dNdDp':
        label = '$\mathrm{d}N\,/\,\mathrm{d}D_{P}$ (nm$^{-1}\,$cm$^{-3}$)'
    elif distType == 'dNdlogDp':
        label = '$\mathrm{d}N\,/\,\mathrm{d}log(D_{P})$ (cm$^{-3}$)'
    elif distType == 'dSdDp':
        label = '$\mathrm{d}S\,/\,\mathrm{d}D_{P}$ (nm$\,$cm$^{-3}$)'
    elif distType == 'dSdlogDp':
        label = '$\mathrm{d}S\,/\,\mathrm{d}log(D_{P})$ (nm$^2\,$cm$^{-3}$)'
    elif distType == 'dVdDp':
        label = '$\mathrm{d}V\,/\,\mathrm{d}D_{P}$ (nm$^2\,$cm$^{-3}$)'
    elif distType == 'dVdlogDp':
        label = '$\mathrm{d}V\,/\,\mathrm{d}log(D_{P})$ (nm$^3\,$cm$^{-3}$)'
    elif distType == 'calibration':
        label = '$\mathrm{d}N\,/\,\mathrm{d}Amp$ (bin$^{-1}\,$cm$^{-3}$)'
    elif distType == 'numberConcentration':
        label = 'Particle number in bin'
    else:
        raise ValueError('%s is not really an option!?!' % distType)
    return label

def get_settings():
    settings = {'wavelength':       {'value': None,
                                     'default': None,
                                     'unit': 'nm'},
                'refractive_index': {'value': None,
                                     'default': None,},
                'mie_result':       {'value': None,
                                     'default': None},
                'particle_density': {'value': 1.8,
                                     'default': 1.8,
                                     'unit': 'g/cc',
                                     'doc': ('Suggested values:\n'
                                             '\tTroposphere: 1.8 g/cc \n'
                                             '\tStratosphere: 2 g/cc')},
                'normalize': {'value': 'None',
                             'default':'None',
                             'unit': 'no unit',
                             'options': ['None', 'STP', 'ambient']},
                'kappa': {'value': None,
                          'default': None,
                          'unit': 'no unit',
                          'doc': ''},
                # 'growth_factor': {'value': None,
                #                   'default':None,
                #                   'unit': 'no unit',
                #                   'doc': 'This should be HygroscopicGrowthFactorDistributions instance.'},
                'growth_distribution': {'value': None,
                                  'default':None,
                                  'unit': 'no unit',
                                  'doc': 'This should be HygroscopicGrowthFactorDistributions instance.'},
                'RH': {'value': None,
                       'default': None,
                       'unit': 'no unit',
                       'doc': 'If None column in housekeeping will be used ... if it exists.'},
                }
    return settings.copy()

def save_netcdf(sizedist, fname, housekeeping = True, value_added_products = True, binunit='nm', tags=[], test=False):
    """Save to netCDF format

    Parameters
    ----------
    fname: str
    housekeeping: bool [True]
        If housekeeping is saved
    value_added_products: bool [True]
        If value added products (currently, norm_fitresults) are going to be saved
    binunit: str
        The units of the diameter of the binedges and centers
    tags: list
    test: bool
        If true the xarray.Dataset is not saved but retured instead.
    """

    sizedist.data.columns.name = 'bincenters'
    sizedist = sizedist.convert2dNdlogDp()
    # distdict = {}
    ds = _xr.Dataset()

    ps = pd.Series(sizedist.bins, name='diameter')
    ps.index.name = 'idx'
    ds['binedges'] = ps

    # atmpy_ps = pd.Series({'type': type(sizedist).__name__})
    # atmpy_ps.index.name = '_atmPy_value'
    ds.attrs['type'] = type(sizedist).__name__
    if type(sizedist).__name__ == 'SizeDist':
        ds['sizedistribution'] = sizedist.data.loc[0, :]
    elif type(sizedist).__name__ == 'SizeDist_TS':
        ds['sizedistribution'] = sizedist.data
        if isinstance(sizedist._data_period, type(None)):
            dp = 'None'
        else:
            dp = sizedist._data_period
        ds.attrs['data_period'] = dp
    elif type(sizedist).__name__ == 'SizeDist_LS':
        sizedist.data.index.name = 'altitude'
        ds['sizedistribution'] = sizedist.data
        df = pd.DataFrame(sizedist.layerbounderies, index=sizedist.layercenters, columns=['low_lim', 'upper_lim'])
        df.index.name = 'altitude'
        ds['layer_boundaries'] = df
    else:
        raise KeyError('not implemented yet... fix it')

    # distdict['_atmPy'] = atmpy_ps

    # ds = _xr.Dataset(distdict)

    ds.attrs['info'] = """Aerosol size distribution data generated with atmPy. For information on variables, e.g. units etc., check the variables attributes."""
    ds.attrs['tags'] = ', '.join(tags)
    ds.binedges.attrs['unit'] = binunit
    ds.bincenters.attrs['unit'] = binunit
    ds.sizedistribution.attrs['moment'] = sizedist.distributionType

    ## Housekeeping
    if housekeeping:
        if sizedist.housekeeping:
            sizedist.housekeeping.data.columns.name = 'hk_columns'
            ds['housekeeping'] = sizedist.housekeeping.data

    ## Normal fit results
    if value_added_products:
        if sizedist._normal_distribution_fits:
            sizedist.normal_distribution_fits.data.columns.name = 'normfit_colums'
            ds['normal_distribution_fits'] = sizedist.normal_distribution_fits.data

    ## Applied corrections
    applied_corrections = []
    if type(sizedist._correct4ambient_LFE_tmp_difference) != type(None):
        applied_corrections.append('flow_in_LFE4temp_difference')

    if type(sizedist._correct4fowrate) != type(None):
        applied_corrections.append('flow_rate')

    ds.attrs['applied_corrections'] = ', '.join(applied_corrections)

    if test:
        return ds
    ds.to_netcdf(fname)
    return

class _Parameter(object):
    def __init__(self, parent, what):
        self._parent = parent
        self._what = what
        self._dict = self._parent._parent._settings[what]

        if 'doc' in self._dict.keys():
            self.__doc__ = self._dict['doc']

    def info(self):
        if 'doc' in self._dict.keys():
            out = self._dict['doc']
        else:
            out = None
        return out

    def __repr__(self):
        return str(self._dict['value'])

    def __bool__(self):
        # try:
        #     bool = self._dict['value'].__bool__()
        # except ValueError:
        return bool(_np.any(self._dict['value']))

    def __len__(self):
        return self._dict['value'].__len__()

    def reset2default(self):
        self._dict['value'] = self._dict['default']


    def _set_value(self, value):
        if 'options' in self._dict.keys():
            if value not in self._dict['options']:
                txt = '{} is not an option for parameter {}. Choose one from {}'.format(value, self.what, self._dict['options'])
                raise ValueError(txt)
        else:
            if type(value).__name__ == '_Parameter':
                value = value._dict['value']
            self._dict['value'] = value

    ###############

    def __add__(self, other):
        return self._dict['value'] + other

    def __radd__(self, other):
        return self._dict['value'] + other

    def __sub__(self, other):
        return self._dict['value'] - other

    def __rsub__(self, other):
        return other - self._dict['value']

    def __mul__(self, other):
        return self._dict['value'] * other

    def __rmul__(self, other):
        return self._dict['value'] * other

    def __truediv__(self, other):
        return self._dict['value'] / other

    def __rtruediv__(self, other):
        return other / self._dict['value']

    ###############
    def copy(self):
        return deepcopy(self)

    @property
    def default_value(self):
        return self._dict['value_default']

    @property
    def description(self):
        return self._dict['description']

    @property
    def value(self):
        return self._dict['value']

# class _Deprecated_SettingOpticalProperty(object):
#     def __init__(self, parent):
#         self._parent = parent
#         self._reset()
#
#     def __repr__(self):
#         out = _all_attributes2string(self)
#         return out
#
#     def _reset(self):
#         self._parent.optical_properties = None
#
#     def _check(self, raise_error = True):
#         dependence_combinations = [['wavelength', 'refractive_index'],  ['mie_results']]
#         passed = []
#         for dependence in dependence_combinations:
#             for dep in dependence:
#                 value = self._parent._settings[dep]['value']
#                 missing = False
#                 if type(value) == type(None):
#                     missing = True
#                 elif not _np.all(~_np.isnan(value)):
#                     missing = True
#
#                 if missing:
#                     # else:
#                     passed.append(False)
#                 else:
#                     passed.append(True)
#         if True in passed:
#             passed = True
#         else:
#             passed = False
#         if not passed:
#             if raise_error:
#                 # txt = 'Parameter {} in optical_property_settings is not set ... do so!'.format(dep)
#                 txt = 'Not enough parameter to calculate optical properties (sizedistribution.optical_properties.parameters).'#.format(dep)
#                 raise ValueError(txt)
#         return passed
#
#     @property
#     def refractive_index(self):
#         return _Parameter(self, 'refractive_index')
#
#     @refractive_index.setter
#     def refractive_index(self, n):
#         if type(n).__name__ in ('int','float'):
#             pass
#         elif type(n).__name__  in ('TimeSeries'):
#             if not _np.array_equal(self._parent.data.index, n.data.index):
#                 n = n.align_to(self)
#             n = n.data
#
#         if type(n).__name__ in ('DataFrame', 'ndarray'):
#             if n.shape[0] != self.data.shape[0]:
#                 txt = """\
# Length of new array has to be the same as that of the size distribution. Use
# sizedistribution.align to align the new array to the appropriate length"""
#                 raise ValueError(txt)
#
#             if not _np.array_equal(self.data.index, n.index):
#                 txt = """\
# The index of the new DataFrame has to be the same as that of the size distribution. Use
# sizedistribution.align to align the index of the new array."""
#                 raise ValueError(txt)
#
#         self._reset()
#         # self._parent._settings['refractive_index']['value'] = n
#         _Parameter(self, 'refractive_index')._set_value(n)
#
#     @property
#     def wavelength(self):
#         # return settings['wavelength']
#         return _Parameter(self, 'wavelength')
#
#     @wavelength.setter
#     def wavelength(self, value):
#         self._reset()
#         # self._parent._settings['wavelength']['value'] = value
#         _Parameter(self, 'wavelength')._set_value(value)



class _Parameters4Reductions(object):
    def __init__(self, parent):

        self._parent = parent
        self._reset_opt_prop()


    def __repr__(self):
        out = _all_attributes2string(self)
        return out

    def _reset_opt_prop(self):
        self._parent._optical_properties = None

    def _reset_hygro(self):
        self._parent._hygroscopicity = None

    def _check_parameter_exists(self, parameters = None, raise_error = True):
        """Check if each parameter in the list of parameters is set to a value other then None.
        If a parameter it self is a list of parameters each list will be treated as a required combination of parameters.
        At least one of the combinations has to be met"""

        passed = False

        for param_list in parameters:
            if not isinstance(param_list, list):
                txt = "This error is due to a recent change. Make sure the each entry is a list of required parameters ... this allows for different combinations of parameters"
                raise ValueError(txt)
            param_option_passed_list = []
            for par in param_list:
                value = self._parent._settings[par]['value']
                if type(value) == type(None):
                    param_option_passed_list.append(False)
                else:
                    param_option_passed_list.append(True)
            if sum(~_np.array(param_option_passed_list)) == 0:
                passed = True

        if not passed:
            txt = "Not enough parameters set. Make sure one of the following parameter combinations are set:\n\t - "
            txt += '\n\t - '.join([', '.join(par) for par in parameters])
            raise ValueError(txt)


        # for param in parameters:
        #     if type(param) != list:
        #         param = [param]
        #         was_list = False
        #     else:
        #         was_list = True
        #
        #     passed_list = []
        #     for par in param:
        #         value = self._parent._settings[par]['value']
        #         exists = True
        #         if type(value) == type(None):
        #             exists = False
        #         # elif type(value).__name__ == 'ndarray':
        #         #     if not _np.all(~_np.isnan(value)):
        #         #         exists = False
        #         # elif hasattr(value, 'data'):
        #         #     if not _np.any(~_np.isnan(value.data)):
        #         #         exists = False
        #         # else:
        #         #     raise ValueError('sorry, programming requried, this type is not allowed')
        #
        #         passed_list.append(exists)
        #
        #     if not _np.any(passed_list):
        #         if raise_error:
        #             if was_list:
        #                 txt = 'One of the Parameters {} or {} needs to be set ... do so!'.format(','.join(param[:-1]), param[-1])
        #             else:
        #                 txt = 'Parameter {} is not set ... do so!'.format(param[0])
        #             raise ValueError(txt)
        #         else:
        #             passed = False

        return passed


        # for par in parameters:
        #     value = self._parent._settings[par]['value']
        #     missing = False
        #     if type(value) == type(None):
        #         missing = True
        #     elif not _np.all(~_np.isnan(value)):
        #         missing = True
        #
        #     if missing:
        #         if raise_error:
        #             txt = 'Parameter {} is not set ... do so!'.format(par)
        #             raise ValueError(txt)
        #         else:
        #             passed = False
        # return passed

    def _check_opt_prop_param_exist(self, raise_error = True):
        return self._check_parameter_exists(parameters= [['wavelength', 'refractive_index'],['mie_result']], raise_error = raise_error)

    def _check_growth_parameters_exist(self, raise_error = True):
        return self._check_parameter_exists(parameters= [['kappa', 'RH'],['growth_distribution', 'RH']], raise_error = raise_error)

    def _check_mixing_ratio_param_exist(self, raise_error = True):
        return self._check_parameter_exists(parameters= [['particle_density']], raise_error = raise_error)

    @property
    def _prop_mie_result(self):
        return _Parameter(self, 'mie_result')

    @_prop_mie_result.setter
    def _prop_mie_result(self,value):
        self._reset_opt_prop()
        _Parameter(self, 'mie_result')._set_value(value)

    @property
    def _prop_refractive_index(self):
        return _Parameter(self, 'refractive_index')

    @_prop_refractive_index.setter
    def _prop_refractive_index(self, n):
        n = align2sizedist(self._parent, n)
        self._reset_opt_prop()
        self._reset_hygro()
        _Parameter(self, 'refractive_index')._set_value(n)

    @property
    def _prop_wavelength(self):
        return _Parameter(self, 'wavelength')

    @_prop_wavelength.setter
    def _prop_wavelength(self, value):
        self._reset_opt_prop()
        _Parameter(self, 'wavelength')._set_value(value)

    @property
    def _prop_kappa(self):
        return _Parameter(self, 'kappa')

    @_prop_kappa.setter
    def _prop_kappa(self, value):
        value = align2sizedist(self._parent, value)
        self._reset_hygro()
        _Parameter(self, 'kappa')._set_value(value)
        _Parameter(self, 'growth_distribution')._set_value(None)

    @property
    def _prop_growth_distribution(self):
        return _Parameter(self, 'growth_distribution')

    @_prop_growth_distribution.setter
    def _prop_growth_distribution(self, value):
        if not isinstance(value, hygroscopicity.HygroscopicGrowthFactorDistributions):
            txt = "Make sure value is of type: atmPy.aerosols.physics.hygroscopicity.HygroscopicGrowthFactorDistributions"
            raise ValueError(txt)
        self._reset_hygro()
        _Parameter(self, 'growth_distribution')._set_value(value)
        _Parameter(self, 'kappa')._set_value(None)

    @property
    def _prop_RH(self):
        return _Parameter(self, 'RH')

    @_prop_RH.setter
    def _prop_RH(self, value):
        self._reset_hygro()
        _Parameter(self, 'RH')._set_value(value)

    @property
    def _prop_particle_density(self):
        return _Parameter(self, 'particle_density')

    @_prop_particle_density.setter
    def _prop_particle_density(self, value):
        """if type timeseries or vertical profile alignment is taken care of"""
        if type(value).__name__ in ['TimeSeries','VerticalProfile']:
            if not _np.array_equal(self._parent.data.index.values, value.data.index.values):
                value = value.align_to(self)
        elif type(value).__name__ not in ['int', 'float']:
            raise ValueError('%s is not an excepted type'%(type(value).__name__))

        self._parent._uptodate_particle_mass_concentration = False
        self._parent._settings['particle_density']['value'] = value

class _Parameters4Reductions_all(_Parameters4Reductions):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        setattr(_Parameters4Reductions_all, 'wavelength', _Parameters4Reductions_all._prop_wavelength)
        setattr(_Parameters4Reductions_all, 'refractive_index', _Parameters4Reductions_all._prop_refractive_index)
        setattr(_Parameters4Reductions_all, 'particle_density', _Parameters4Reductions_all._prop_particle_density)
        setattr(_Parameters4Reductions_all, 'kappa', _Parameters4Reductions_all._prop_kappa)
        setattr(_Parameters4Reductions_all, 'growth_distribution', _Parameters4Reductions_all._prop_growth_distribution)
        setattr(_Parameters4Reductions_all, 'RH', _Parameters4Reductions_all._prop_RH)
        setattr(_Parameters4Reductions_all, 'mie_result', _Parameters4Reductions_all._prop_mie_result)

class _Parameters4Reductions_opt_prop(_Parameters4Reductions):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        setattr(_Parameters4Reductions_opt_prop, 'wavelength', _Parameters4Reductions_opt_prop._prop_wavelength)
        setattr(_Parameters4Reductions_opt_prop, 'refractive_index', _Parameters4Reductions_opt_prop._prop_refractive_index)
        setattr(_Parameters4Reductions_opt_prop, 'mie_result', _Parameters4Reductions_opt_prop._prop_mie_result)


class _Parameters4Reductions_hygro_growth(_Parameters4Reductions):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        setattr(_Parameters4Reductions_hygro_growth, 'refractive_index',    _Parameters4Reductions_hygro_growth._prop_refractive_index)
        setattr(_Parameters4Reductions_hygro_growth, 'particle_density',    _Parameters4Reductions_hygro_growth._prop_particle_density)
        setattr(_Parameters4Reductions_hygro_growth, 'kappa',               _Parameters4Reductions_hygro_growth._prop_kappa)
        setattr(_Parameters4Reductions_hygro_growth, 'growth_distribution', _Parameters4Reductions_hygro_growth._prop_growth_distribution)
        setattr(_Parameters4Reductions_hygro_growth, 'RH',                  _Parameters4Reductions_hygro_growth._prop_RH)


# class _SettingHygroscopicGrowth(object):
#     def __init__(self, parent):
#         self._parent = parent
#         # self._reset()
#
#     def _check(self, raise_error = True):
#         dependence = ['refractive_index']
#         passed = True
#         for dep in dependence:
#             value = self._parent._settings[dep]['value']
#
#             missing = False
#             if type(value) == type(None):
#                 missing = True
#             elif not _np.all(~_np.isnan(value)):
#                 missing = True
#
#             if missing:
#                 if raise_error:
#                     txt = 'Parameter {} in hygroscopic_growth_settings is not set ... do so!'.format(dep)
#                     raise ValueError(txt)
#                 else:
#                     passed = False
#         return passed
#
#     @property
#     def _prop_refractive_index(self):
#         return _Parameter(self, 'refractive_index')
#
#     @_prop_refractive_index.setter
#     def _prop_refractive_index(self, n):
#         if type(n).__name__ in ('int', 'float'):
#             pass
#         elif type(n).__name__ in ('TimeSeries'):
#             if not _np.array_equal(self._parent.data.index, n.data.index):
#                 n = n.align_to(self)
#             n = n.data
#
#         if type(n).__name__ in ('DataFrame', 'ndarray'):
#             if n.shape[0] != self._parent.data.shape[0]:
#                 txt = """\
# Length of new array has to be the same as that of the size distribution. Use
# sizedistribution.align to align the new array to the appropriate length"""
#                 raise ValueError(txt)
#
#             if not _np.array_equal(self._parent.data.index, n.index):
#                 txt = """\
# The index of the new DataFrame has to the same as that of the size distribution. Use
# sizedistribution.align to align the index of the new array."""
#                 raise ValueError(txt)
#
#         # self._reset()
#         # self._parent._settings['refractive_index']['value'] = n
#         _Parameter(self, 'refractive_index')._set_value(n)
#
def _all_attributes2string(obj):
    att_list = []
    max_len = 0
    for i in dir(obj):
        if i[0] == '_':
            continue
        if len(i) > max_len:
            max_len = len(i)

    for i in dir(obj):
        if i[0] == '_':
            continue
        att_list.append('{i:<{max_len}}:  {value}'.format(i=i, max_len=max_len + 1, value=getattr(obj, i)))
    out = '\n'.join(att_list)
    return out

class _Properties(object):
    def __init__(self, parent):
        self._parent = parent

    def __repr__(self):
        out = _all_attributes2string(self)
        return out

    @property
    def particle_density(self):
        return _Parameter(self, 'particle_density')

    @particle_density.setter
    def particle_density(self, value):
        """if type timeseries or vertical profile alignment is taken care of"""
        if type(value).__name__ in ['TimeSeries','VerticalProfile']:
            if not _np.array_equal(self._parent.data.index.values, value.data.index.values):
                value = value.align_to(self)
        elif type(value).__name__ not in ['int', 'float']:
            raise ValueError('%s is not an excepted type'%(type(value).__name__))

        self._parent._uptodate_particle_mass_concentration = False
        self._parent._settings['particle_density']['value'] = value


class _NormalDistributionFitRes_VP(_vertical_profile.VerticalProfile):
    def __init__(self, parent_dist, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent_dist = parent_dist

    def plot_fitres(self, show_width = True, show_amplitude = True, ax = None, **kwargs):
        """ Plots the results from fit_normal

        Arguments
        ---------
        amp: bool.
            if the amplitude is to be plotted
        """
        if ax == None:
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()

        if show_width:
            a.fill_betweenx(self.data.index, self.data.Sigma_high, self.data.Sigma_low,
                           color=_colors[0],
                           alpha=0.5,
                           )

        a.plot(self.data.Pos, self.data.index, **kwargs)
        g = a.get_lines()[-1]
        g.set_label('Center')
        a.legend(loc=2)

        a.set_xlabel('Particle diameter (nm)')
        a.set_ylabel('Altitude (m)')


        if show_amplitude:
            a2 = a.twiny()
            a2.plot(self.data.Amp, self.data.index)
            # self.data_fit_normal.Amp.plot(ax=a2, color=plt_tools.color_cycle[1], linewidth=2)
            g = a2.get_lines()[-1]
            g.set_color(_colors[1])
            g.set_label('Amplitude of norm. dist.')
            a2.legend()
            a2.set_ylabel('Amplitude - %s' % (get_label(self.parent_dist.distributionType)))
            a2.set_xscale('log')
        else:
            a2 = None
        return f, a, a2



class _NormalDistributionFitRes_TS(_timeseries.TimeSeries):
    def __init__(self, parent_dist, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.parent_dist = parent_dist

    def plot_fitres(self, show_width = True, show_amplitude = True, ax = None, **kwargs):
        """ Plots the results from fit_normal"""

        if ax:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()
        data = self.data.dropna()
        if show_width:
            a.fill_between(data.index, data.Sigma_high, data.Sigma_low,
                           color=plt_tools.color_cycle[0],
                           alpha=0.5,
                           )
        a.plot(data.index.values, data.Pos.values, **kwargs)
        # data.Pos.plot(ax=a, color=plt_tools.color_cycle[0], linewidth=2, label='center')
        a.legend(loc=2)
        a.set_ylabel('Particle diameter (nm)')
        a.set_xlabel('')

        if show_amplitude:
            a2 = a.twinx()
            # data.Amp.plot(ax=a2, color=plt_tools.color_cycle[1], linewidth=2, label='amplitude')
            a2.plot(data.index.values, data.Amp.values, color=plt_tools.color_cycle[1], linewidth=2, label='amplitude')
            a2.legend()
            a2.set_ylabel('Amplitude - %s' % (get_label(self.parent_dist.distributionType)))
            a2.set_yscale('log')
        else:
            a2 = None
        f.autofmt_xdate()
        return f, a, a2



class SizeDist(object):
    """
    Object defining a log normal aerosol size distribution


    Arguments
    ----------
    bincenters:         NumPy array, optional
                        this is if you actually want to pass the bincenters, if False they will be calculated
    distributionType:
                        log normal: 'dNdlogDp','dSdlogDp','dVdlogDp'
                        natural: 'dNdDp','dSdDp','dVdDp'
                        number: 'dNdlogDp', 'dNdDp', 'numberConcentration'
                        surface: 'dSdlogDp','dSdDp'
                        volume: 'dVdlogDp','dVdDp'
    data:   pandas dataFrame, optional
            None, will generate an empty pandas data frame with columns defined by bins
             - pandas dataFrame with
                 - column names (each name is something like this: '150-200')
                 - index is time (at some point this should be arbitrary, convertable to altitude for example?)
    unit conventions:
             - diameters: nanometers
             - flowrates: cc (otherwise, axis label need to be adjusted an caution needs to be taken when dealing is AOD)



    Notes
    ------
    * Diameters are specified in nanometers

    """
    # todo: write setters and getters for bins and bincenter, so when one is changed the otherone is automatically
    #  changed too
    def __init__(self, data, bins, distType,
                 # bincenters=False,
                 # fixGaps=False
                 ):

        if type(data).__name__ == 'NoneType':
            self._data = pd.DataFrame()
        else:
            self._data = data

        self._settings = get_settings()
        self._optical_properties = None
        self._hygroscopicity = None

        # self.optical_properties_settings = _SettingOpticalProperty(self)
        self.parameters4reductions = _Parameters4Reductions_all(self)
        # self.hygroscopic_growth_settings = _SettingHygroscopicGrowth(self)
        # self.properties = _Properties(self)

        self.bins = bins
        # self.__index_of_refraction = None
        # self._growth_factor = None
        self.__particle_number_concentration = None
        self.__particle_mass_concentration = None
        self.__particle_surface_concentration = None
        self.__particle_volume_concentration = None
        self._submicron_volume_ratio = None
        self._housekeeping = None

        self.distributionType = distType
        self._sup_opt_wl = None #todo: deprecated?!?
        self.__sup_opt_aod = False #todo: deprecated?!?
        self.particle_number_concentration_outside_range = None
        self._update()
        self._is_reduced_to_pt = False
        self._mode_analysis = None
        self._normal_distribution_fits = None

        self._correct4ambient_LFE_tmp_difference = None
        self._correct4fowrate = None
        self._mean_diameter = None


        # if fixGaps:
        #     self.fillGaps()

    def __mul__(self, other):
        dist  = self.copy()
        dist.data *= other
        if dist.particle_number_concentration_outside_range:
            dist.particle_number_concentration_outside_range *= other

        dist._update()
        return dist

    def __truediv__(self, other):
        self.data /= other
        self.particle_number_concentration_outside_range /= other
        self._update()
        return self

    def __add__(self, other):
        if type(other) == type(self):
            if self.distributionType != other.distributionType:
                raise ValueError('size distributions have to have the same distributionType')
            if self.data.shape != other.data.shape:
                raise ValueError('size distributions have to have the same shape')
            self.data = self.data + other.data
            # self.particle_number_concentration_outside_range = self.particle_number_concentration_outside_range + other.particle_number_concentration_outside_range
        else:
            self.data += other
            self.particle_number_concentration_outside_range += other
        self._update()
        return self

    def __sub__(self, other):
        self.data -= other
        self.particle_number_concentration_outside_range -= other
        self._update()
        return self


    # mode_analysis = modes.ModeAnalysis

    @property
    def data(self):
        return self._data

    @data.setter
    def data(self, value):
        self._data = value
        self._update()

    @property
    def mode_analysis(self):
        if not self._mode_analysis:
            self._mode_analysis = modes.ModeAnalysis(self)
        return self._mode_analysis

    @property
    def DEPRECATEDoptical_properties(self):
        if not self._optical_properties:
            self._optical_properties = optical_properties.size_dist2optical_properties(self, aod = self.__sup_opt_aod, noOfAngles=100)
        return self._optical_properties

    @property
    def optical_properties(self):
        if not self._optical_properties:
            self._optical_properties = optical_properties.OpticalProperties(self)
        return self._optical_properties

    @property
    def hygroscopicity(self):
        if not self._hygroscopicity:
            self._hygroscopicity = hygroscopicity.HygroscopicityAndSizedistributions(self)
        return self._hygroscopicity

    # @optical_properties.setter
    # def optical_properties(self, value):
    #     self.__optical_properties = value

    # @property
    # def physical_property_density(self):
    #     """setter: if type _timeseries or vertical profile alignment is taken care of"""
    #     return self.__physical_property_density
    #
    # @physical_property_density.setter
    # def physical_property_density(self, value):
    #     """if type timeseries or vertical profile alignment is taken care of"""
    #     if type(value).__name__ in ['TimeSeries','VerticalProfile']:
    #         if not _np.array_equal(self.data.index.values, value.data.index.values):
    #             value = value.align_to(self)
    #     elif type(value).__name__ not in ['int', 'float']:
    #         raise ValueError('%s is not an excepted type'%(type(value).__name__))
    #     self.__physical_property_density = value

    @property
    def housekeeping(self):
        return self._housekeeping

    @housekeeping.setter
    def housekeeping(self, value):
        if value:
            self._housekeeping = value.align_to(self)
        else:
            self._housekeeping = value

    @property
    def bins(self):
        return self.__bins

    @bins.setter
    def bins(self,array):
        # bins_st = array.astype(int).astype(str)
        # col_names = []
        self.__bins = array
        self.__bincenters = (array[1:] + array[:-1]) / 2.
        self.__binwidth = (array[1:] - array[:-1])
        # self.data.columns = _np.round(self.bincenters, 0).astype(_np.float32)
        self.data.columns = self.bincenters
        self.data.columns.name = 'bincenters'

    @property
    def bincenters(self):
        return self.__bincenters

    @property
    def binwidth(self):
        return self.__binwidth
    
    @property
    def mean_diameter_geometric(self):
        if isinstance(self._mean_diameter, type(None)):
            md = self.convert2numberconcentration().data.apply(lambda row:(row.index ** row.astype(_np.float128)).prod()**(1/row.sum()), axis = 1)
            self._mean_diameter = md
        return self._mean_diameter
    

    @property
    def normal_distribution_fits(self):
        if type(self._normal_distribution_fits) == type(None):
            self._normal_distribution_fits = fit_normal_distribution2sizedist(self)
        return self._normal_distribution_fits

    @normal_distribution_fits.setter
    def normal_distribution_fits(self,value):
        self._normal_distribution_fits = value

    @property
    def sup_optical_properties_wavelength(self):
        return self._sup_opt_wl

    @sup_optical_properties_wavelength.setter
    def sup_optical_properties_wavelength(self, data):
        self._optical_properties = None
        self._sup_opt_wl = data

#     @property
#     def index_of_refraction(self):
#         """In case of setting the value and value is TimeSeries it will be aligned to the time series of the
#         size distribution.
#         """
#         return self.__index_of_refraction
#
#     @index_of_refraction.setter
#     def index_of_refraction(self,n):
#         self.__optical_properties = None
#         if type(n).__name__ in ('int','float'):
#             pass
#         elif type(n).__name__  in ('TimeSeries'):
#             if not _np.array_equal(self.data.index, n.data.index):
#                 n = n.align_to(self)
#             n = n.data
#
#         if type(n).__name__ in ('DataFrame', 'ndarray'):
#             if n.shape[0] != self.data.shape[0]:
#                 txt = """\
# Length of new array has to be the same as that of the size distribution. Use
# sizedistribution.align to align the new array to the appropriate length"""
#                 raise ValueError(txt)
#
#             if not _np.array_equal(self.data.index, n.index):
#                 txt = """\
# The index of the new DataFrame has to the same as that of the size distribution. Use
# sizedistribution.align to align the index of the new array."""
#                 raise ValueError(txt)
#
#         self.__index_of_refraction = n
        # elif self.__index_of_refraction:
        #     txt = """Security stop. This is to prevent you from unintentionally changing this value.
        #     The index of refraction is already set to %.2f, either by you or by another function, e.g. apply_hygro_growth.
        #     If you really want to change the value do it by setting the __index_of_refraction attribute."""%self.index_of_refraction
        #     raise ValueError(txt)

    # @property
    # def growth_factor(self):
    #     return self._growth_factor

    @property
    def particle_number_concentration(self):
        if not _np.any(self.__particle_number_concentration) or not self._uptodate_particle_number_concentration:
            self.__particle_number_concentration = self._get_particle_concentration()
        return self.__particle_number_concentration

    @property
    def particle_mean_diameter(self):
        if not _np.any(self._particle_mean_diameter):
            self._particle_mean_diameter = self._get_particle_mean_diameter()
        return self._particle_mean_diameter

    def _get_particle_mean_diameter(self):
        def mean_d(row):
            if row.sum() == 0:
                return _np.nan
            return (row.index.values * row).sum() / row.sum()
        md = self.data.apply(mean_d, axis=1)
        df = pd.DataFrame(md, columns=['mean d (nm)'])
        return df

    @property
    def particle_mass_concentration(self):
        if not _np.any(self.__particle_mass_concentration) or not self._uptodate_particle_mass_concentration:
            self.__particle_mass_concentration = self._get_mass_concentration()
            self._uptodate_particle_mass_concentration = True
        return self.__particle_mass_concentration


    @property
    def particle_mass_mixing_ratio(self):
        raise AttributeError('sorry does not exist for a single size distribution')

    @property
    def particle_number_mixing_ratio(self):
        raise AttributeError('sorry does not exist for a single size distribution')

    @property
    def particle_volume_concentration(self):
        if not _np.any(self.__particle_volume_concentration) or not self._uptodate_particle_volume_concentration:
            self.__particle_volume_concentration = self._get_volume_concentration()
            self._uptodate_particle_volume_concentration = True
        return self.__particle_volume_concentration

    @property
    def particle_surface_concentration(self):
        if not _np.any(self.__particle_surface_concentration) or not self._uptodate_particle_surface_concentration:
            self.__particle_surface_concentration = self._get_surface_concentration()
            self._uptodate_particle_surface_concentration = True
        return self.__particle_surface_concentration

    @property
    def submicron_volume_ratio(self):
        if not self._submicron_volume_ratio:
            self._submicron_volume_ratio = self.zoom_diameter(end=1000).particle_volume_concentration / self.particle_volume_concentration
        return self._submicron_volume_ratio

    re_bin = _diameter_binning.re_bin

    def reduce2temp_press_ambient(self, tmp_is = 'auto', tmp_is_column = 'Temperature_instrument', press_is_column = 'Pressure_Pa'):
        """This function corrects the particles concentrations to ambient conditions. This is necessary if the
        temperature of the instrument is different then ambient. When the instrument is adjusting the flow to a constant
        rate it will be at the instrument temperature not ambient -> correction required
        tmp in C
        press in hPa"""
        dist = self.copy()
        if dist._is_reduced_to_pt:
            raise TypeError('Already is reduced to ambient conditions')
        data = dist.data
        press_is = dist.housekeeping.data[press_is_column]
        if tmp_is == 'auto':
            tmp_is = dist.housekeeping.data[tmp_is_column] + 273.15
        else:
            tmp_is+=273.15

        press_shall = dist.housekeeping.data['Pressure_Pa']
        tmp_shall = dist.housekeeping.data['Temperature'] + 273.15

        newdata = data.mul(tmp_is, axis=0).truediv(press_is, axis=0).truediv(tmp_shall, axis=0).mul(press_shall, axis=0)
        dist.data = newdata
        dist._update()
        dist._is_reduced_to_pt = True
        return dist

    def reduce2temp_press_standard(self, tmp_is = 'auto', tmp_is_column = 'Temperature_instrument', press_is_column = 'Pressure_Pa'):
        """tmp in C
        press in hPa"""
        dist = self.copy()
        if dist._is_reduced_to_pt:
            raise TypeError('Already is reduced to ambient conditions')
        data = dist.data
        press_is = dist.housekeeping.data[press_is_column]

        if tmp_is == 'auto':
            tmp_is = dist.housekeeping.data[tmp_is_column] + 273.15
        else:
            tmp_is+=273.15

        press_shall = 1000
        tmp_shall = 273.15

        newdata = data.mul(tmp_is, axis=0).truediv(press_is, axis=0).truediv(tmp_shall, axis=0).mul(press_shall, axis=0)
        dist.data = newdata
        dist._update()
        dist._is_reduced_to_pt = True
        return dist
    

    def correct4flowrate(self, flowrate):
        """This simply normalizes to to provided flow rate. In principle this could be the flow rate reported in the
        the housekeeping file ... and in other instruments then POPS it probably shoud... However in POPS this value is
        questionable and the set_point is a more usefull value"""

        sizedist = self.copy()
        if type(sizedist._correct4fowrate) == type(None):
            sizedist *= 1/flowrate
            sizedist._correct4fowrate = flowrate
        else:
            txt = "You can't apply this correction twice!"
            raise ValueError(txt)
            # _warnings.warn(txt)
        return sizedist

    def correct4ambient_LFE_tmp_difference(self):
        """corrects for temperature differences between ambient and instrument.
        The Problem is that the instrument samples at a constant flow at the temperature of
        the laminar flow element not the ambient temperatrue ... this corrects for it
        Make sure your housekeeping has a column named Temperature and one named Temperature_instrument"""
        sizedist = self.copy()
        if type(sizedist._correct4ambient_LFE_tmp_difference) == type(None):
            tcorr = (sizedist.housekeeping.data.Temperature_instrument + 273.15) / (sizedist.housekeeping.data.Temperature + 273.15)
            sizedist.data = sizedist.data.apply(lambda x: x * tcorr, axis=0)
            sizedist._update()
            sizedist._correct4ambient_LFE_tmp_difference = tcorr
        else:
            txt = "You can't apply this correction twice!"
            raise ValueError(txt)
            # _warnings.warn(txt)
        return sizedist

    # def deprecated_apply_hygro_growth(self, kappa, RH, how ='shift_bins', adjust_refractive_index = True):
    #     """Note kappa values are !!NOT!! aligned to self in case its timesersies
    #     how: string ['shift_bins', 'shift_data']
    #         If the shift_bins the growth factor has to be the same for all lines in
    #         data (important for timeseries and vertical profile.
    #         If gf changes (as probably the case in TS and LS) you want to use
    #         'shift_data'
    #     """
    #
    #     self.parameters4reductions._check_growth_parameters_exist()
    #
    #
    #     dist_g = self.convert2numberconcentration()
    #     # dist_g = dist_g.copy()
    #     # pdb.set_trace()
    #
    #     if type(kappa).__name__ == 'TimeSeries':
    #         kappa = kappa.data.iloc[:, 0].values
    #     # gf,n_mix = hg.kappa_simple(kappa, RH, refractive_index= dist_g.index_of_refraction)
    #     gf,n_mix = hygroscopicity.kappa_simple(kappa, RH, refractive_index= dist_g.parameters4reductions.refractive_index)
    #     # pdb.set_trace()
    #
    #     if how == 'shift_bins':
    #         if not isinstance(gf, (float,int)):
    #             txt = '''If how is equal to 'shift_bins' RH has to be of type int or float.
    #             It is %s'''%(type(RH).__name__)
    #             raise TypeError(txt)
    #     if type(self).__name__ == 'SizeDist_LS':
    #         if type(gf) in (float, int):
    #             nda = _np.zeros(self.data.index.shape)
    #             nda[:] = gf
    #             gf = nda
    #         gf = _vertical_profile.VerticalProfile(pd.DataFrame(gf, index = self.data.index))
    #
    #     elif type(self).__name__ == 'SizeDist_TS':
    #         if type(gf).__name__ in ('float', 'int', 'float64'):
    #             nda = _np.zeros(self.data.index.shape)
    #             nda[:] = gf
    #             gf = nda
    #         # import pdb; pdb.set_trace()
    #         gf = _timeseries.TimeSeries(pd.DataFrame(gf, index = self.data.index))
    #         gf._data_period = self._data_period
    #     # pdb.set_trace()
    #     dist_g = dist_g.deprecated_apply_growth(gf, how = how)
    #
    #     # pdb.set_trace()
    #     if how == 'shift_bins':
    #         dist_g.parameters4reductions.refractive_index = n_mix
    #     elif how == 'shift_data':
    #         if adjust_refractive_index:
    #             if type(n_mix).__name__ in ('float', 'int', 'float64'):
    #                 nda = _np.zeros(self.data.index.shape)
    #                 nda[:] = n_mix
    #                 n_mix = nda
    #             # import pdb;
    #             # pdb.set_trace()
    #             df = pd.DataFrame(n_mix)
    #             df.columns = ['index_of_refraction']
    #             df.index = dist_g.data.index
    #             # pdb.set_trace()
    #             dist_g.parameters4reductions.refractive_index = df
    #         else:
    #             # print('no')
    #             dist_g.parameters4reductions.refractive_index = self.parameters4reductions.refractive_index
    #
    #     # pdb.set_trace()
    #     # dist_g._growth_factor = pd.DataFrame(gfsdf, index = dist_g.data.index, columns = ['Growth_factor'])
    #     # pdb.set_trace()
    #     return dist_g
    #
    #
    # def deprecated_apply_growth(self, growth_factor, how ='auto'):
    #     """Note this does not adjust the refractive index according to the dilution!!!!!!!"""
    #     # pdb.set_trace()
    #     if how == 'auto':
    #         if isinstance(growth_factor, (float, int)):
    #             how = 'shift_bins'
    #         else:
    #             how = 'shift_data'
    #     # pdb.set_trace()
    #     dist_g = self.convert2numberconcentration()
    #     # pdb.set_trace()
    #     if how == 'shift_bins':
    #         if not isinstance(growth_factor, (float, int)):
    #             txt = '''If how is equal to 'shift_bins' the growth factor has to be of type int or float.
    #             It is %s'''%(type(growth_factor).__name__)
    #             raise TypeError(txt)
    #         dist_g.bins = dist_g.bins * growth_factor
    #
    #     elif how == 'shift_data':
    #         if isinstance(growth_factor, (float, int)):
    #             pass
    #         # elif type(growth_factor).__name__ == 'ndarray':
    #         #     growth_factor = _timeseries.TimeSeries(growth_factor)
    #
    #         # elif type(growth_factor).__name__ == 'Series':
    #         #     growth_factor = _timeseries.TimeSeries(pd.DataFrame(growth_factor))
    #         #
    #         # elif type(growth_factor).__name__ == 'DataFrame':
    #         #     growth_factor = _timeseries.TimeSeries(growth_factor)
    #
    #         elif type(growth_factor).__name__ == 'VerticalProfile':
    #             pass
    #
    #         elif type(growth_factor).__name__ == 'TimeSeries':
    #             if growth_factor._data_period == None:
    #                 growth_factor._data_period = self._data_period
    #             growth_factor = growth_factor.align_to(dist_g)
    #         else:
    #             txt = 'Make sure type of growthfactor is int,float,TimeSeries, or Series. It currently is: %s.'%(type(growth_factor).__name__)
    #             raise TypeError(txt)
    #         try:
    #             growth_max = float(_np.nanmax(growth_factor.data))
    #         except AttributeError:
    #             growth_max = growth_factor
    #
    #         test = dist_g._hygro_growht_shift_data(dist_g.data.values[0], dist_g.bins, growth_max, ignore_data_nan = True)
    #         bin_num = test['data'].shape[0]
    #         data_new = _np.zeros((dist_g.data.shape[0],bin_num), dtype= object) # this has to be of type object, so the sum is nan when all are nan, otherwise it would be 0
    #         data_new[:] = _np.nan
    #         #todo: it would be nicer to have _hygro_growht_shift_data take the TimeSeries directly
    #         gf = growth_factor.data.values.transpose()[0]
    #         for e,i in enumerate(dist_g.data.values):
    #             out = dist_g._hygro_growht_shift_data(i, dist_g.bins, gf[e])
    #             dt = out['data']
    #             diff = bin_num - dt.shape[0]
    #             diff_nans = _np.zeros(diff, dtype= object)
    #             diff_nans[:] = _np.nan
    #             dt = _np.append(dt, diff_nans)
    #             data_new[e] = dt
    #         df = pd.DataFrame(data_new)
    #         df.index = dist_g.data.index
    #
    #         # if type(dist_g)
    #
    #         # pdb.set_trace()
    #         dist_g_new = SizeDist(df, test['bins'], dist_g.distributionType)
    #         # pdb.set_trace()
    #         dist_g_new.parameters4reductions.refractive_index = dist_g.parameters4reductions.refractive_index
    #         # pdb.set_trace()
    #         dist_g_new.parameters4reductions.wavelength =       dist_g.parameters4reductions.wavelength
    #         # pdb.set_trace()
    #         # dist_g_new.optical_properties_settings.wavelength =       dist_g.optical_properties_settings.wavelength
    #         # if type(growth_factor).__name__ == 'TimeSeries':
    #         #     dp = dist_g._data_period
    #         #     dist_g._data_period = dp
    #
    #     else:
    #         txt = '''How has to be either 'shift_bins' or 'shift_data'.'''
    #         raise ValueError(txt)
    #
    #     return dist_g_new


    # def grow_particles(self, shift=1):
    #     """This function shifts the data by "shift" columns to the right
    #     Argurments
    #     ----------
    #     shift: int.
    #         number of columns to shift.
    #
    #     Returns
    #     -------
    #     New dist_LS instance
    #     Growth ratio (mean,std) """
    #
    #     dist_grow = self.copy()
    #     gf = dist_grow.bincenters[shift:] / dist_grow.bincenters[:-shift]
    #     gf_mean = gf.mean()
    #     gf_std = gf.std()
    #
    #     shape = dist_grow.data.shape[1]
    #     dist_grow.data[:] = 0
    #     dist_grow.data.iloc[:, shift:] = self.data.values[:, :shape - shift]
    #
    #     return dist_grow, (gf_mean, gf_std)

    # def calculate_optical_properties(self, wavelength, n):
    #     out = optical_properties.size_dist2optical_properties(self, wavelength, n)
    #     return out

    # todo: this function appears multiple times, can easily be inherited
    # def calculate_optical_properties(self, wavelength, n = None, AOD = False, noOfAngles=100):
    #     if not _np.any(n):
    #         n = self.index_of_refraction
    #     if not _np.any(n):
    #         txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #         raise ValueError(txt)
    #     out = optical_properties.size_dist2optical_properties(self, wavelength, n, aod = AOD, noOfAngles=noOfAngles)
    #     opt_properties = optical_properties.OpticalProperties(out, parent = self)
    #     # opt_properties.wavelength = wavelength #should be set in OpticalProperty class
    #     # opt_properties.index_of_refractio = n
    #     #opt_properties.angular_scatt_func = out['angular_scatt_func']  # This is the formaer phase_fct, but since it is the angular scattering intensity, i changed the name
    #     # opt_properties.parent_dist = self
    #     return opt_properties

    def fillGaps(self, scale=1.1):
        """Note: This function is purly esteticall and should be removed since it can also create errors ...
        Finds gaps in dataset (e.g. when instrument was shut of) and fills them with zeros.

        It adds one line of zeros to the beginning and one to the end of the gap. 
        Therefore the gap is visible as zeros instead of the interpolated values

        Parameters
        ----------
        scale:  float, optional
                This is a scale.
        """
        _warnings.warn('This function is deprecated use close_gaps instead ... unless you like thisone better')
        diff = self.data.index[1:].values - self.data.index[0:-1].values
        threshold = _np.median(diff) * scale
        where = _np.where(diff > threshold)[0]
        if len(where) != 0:
            _warnings.warn('The dataset provided had %s gaps' % len(where))
            gap_start = self.data.index[where]
            gap_end = self.data.index[where + 1]
            for gap_s in gap_start:
                self.data.loc[gap_s + threshold] = _np.zeros(self.bincenters.shape)
            for gap_e in gap_end:
                self.data.loc[gap_e - threshold] = _np.zeros(self.bincenters.shape)
            self.data = self.data.sort_index()
        return


    def plot(self,
             showMinorTickLabels=True,
             removeTickLabels=["700", "900"],
             fit_res=True,
             fit_res_scale = 'log',
             ax=None,
             **kwargs
             ):
        """
        Plots and returns f,a (figure, axis).

        Arguments
        ---------
        showMinorTickLabels: bool [True], optional
            if minor tick labels are labled
        removeTickLabels: list of string ["700", "900"], optional
            list of tick labels aught to be removed (in case there are overlapping)
        fit_res: bool [True], optional
            allows plotting of fitresults if fit_normal was previously executed
        fit_res: string
            If fit_normal was done using log = False, you want to set this to linear!
        ax: axis object [None], optional
            option to provide axis to plot on

        Returns
        -------
        Handles to the figure and axes of the figure.


        """
        if type(ax).__name__ in _axes_types:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()

        g, = a.plot(self.bincenters, self.data.iloc[0,:], **kwargs)
        g.set_drawstyle('steps-mid')

        a.set_xlabel('Particle diameter (nm)')

        label = get_label(self.distributionType)
        a.set_ylabel(label)
        a.set_xscale('log')

        if fit_res:
            if 'data_fit_normal' in dir(self):
                amp, pos, sigma = self.data_fit_normal.values[0, :3]
                if fit_res_scale == 'log':
                    normal_dist = math_functions.gauss(_np.log10(self.bincenters), amp, _np.log10(pos), sigma)
                elif fit_res_scale =='linear':
                    normal_dist = math_functions.gauss(self.bincenters, amp, pos, sigma)
                else:
                    txt = '"fit_res_scale has to be either log or linear'
                    raise ValueError(txt)
                a.plot(self.bincenters, normal_dist,
                       label='fit with norm. dist.')
                a.legend()

        return f, a


    def convert2dNdDp(self):
        return self._convert2otherDistribution('dNdDp')

    def convert2dNdlogDp(self):
        return self._convert2otherDistribution('dNdlogDp')

    def convert2dSdDp(self):
        return self._convert2otherDistribution('dSdDp')

    def convert2dSdlogDp(self):
        return self._convert2otherDistribution('dSdlogDp')

    def convert2dVdDp(self):
        return self._convert2otherDistribution('dVdDp')

    def convert2dVdlogDp(self):
        return self._convert2otherDistribution('dVdlogDp')

    def convert2numberconcentration(self):
        return self._convert2otherDistribution('numberConcentration')

    def copy(self):
        return deepcopy(self)

    def extend_bin_range(self, newlimit):
        """Extends the bin range. This will only work if bins are log-equally spaced. Only tested for SizeDist (not for TS or LS). Currently works on the small
        diameter side only ... easy programming will allow for larger diameters too

        Parameters
        ----------
        newlimit: float
            The new lower diameter limit.

        Returns
        -------
        size distribution instance with extra bins. Data filled with nan."""

        self = self.copy()
        newmin = newlimit
        lbins = _np.log10(self.bins)
        lbinsdist = _np.unique(lbins[1:] - lbins[:-1])
        if lbinsdist.shape[0] != 1:
            raise ValueError(
                'The binwidth should be identical thus this value should be one ... programming needed to deal with exceptions')
        while lbins[0] > _np.log10(newmin):
            lbins = _np.append(_np.array([lbins[0] - lbinsdist]), lbins)
            newbins = 10 ** lbins
            newcenter = (newbins[1] + newbins[0]) / 2.
            self.data[newcenter] = _np.nan

        self.data = self.data.transpose().sort_index().transpose()

        self.bins = newbins
        return self

    def extrapolate_size_dist(self, newlimit):
        """Extrapolates the size distribution range assuming a nomal distributed aerosol mode. This will only work if bins are log-equally spaced. Only tested for SizeDist (not for TS or LS). Currently works on the small
        diameter side only ... easy programming will allow for larger diameters too

        Parameters
        ----------
        newlimit: float
            The new lower diameter limit.

        Returns
        -------
        size distribution instance with extra bins. Data filled with results from fitting with normal distribution."""

        sizedist = self.copy()
        sizedist._update()
        sizedist = sizedist.extend_bin_range(newlimit)
        amp, pos, sigma = sizedist.data_fit_normal.values[0, :3]
        normal_dist = math_functions.gauss(_np.log10(sizedist.bincenters), amp, _np.log10(pos), sigma)
        sizedist.data.loc[0, :][_np.isnan(sizedist.data.loc[0, :])] = normal_dist[_np.isnan(sizedist.data.loc[0, :])]
        return sizedist

    def grow_sizedistribution(self, growthfactor):
        return hygroscopicity.apply_growth2sizedist(self, growthfactor)

    def save_csv(self, fname, header=True):
        if header:
            raus = open(fname, 'w')
            raus.write('bins = %s\n' % self.bins.tolist())
            raus.write('distributionType = %s\n' % self.distributionType)
            raus.write('objectType = %s\n' % (type(self).__name__))
            raus.write('#\n')
            raus.close()
        self.data.to_csv(fname, mode='a')
        return

    save_netcdf = save_netcdf

    def save_hdf(self, hdf, variable_name = None, info = '', force = False):

        if variable_name:
            table_name = '/atmPy/aerosols/sizedistribution/'+variable_name
            if table_name in hdf.keys():
                    if not force:
                        txt = 'Table name (variable_name) exists. If you want to overwrite it set force to True.'
                        raise KeyError(txt)
        else:
            e = 0
            while 1:
                table_name = '/atmPy/aerosols/sizedistribution/'+ type(self).__name__ + '_%.3i'%e
                if table_name in hdf.keys():
                    e+=1
                else:
                    break


        hdf.put(table_name, self.data)

        storer = hdf.get_storer(table_name)

        attrs = {}
        attrs['variable_name'] = variable_name
        attrs['info'] = info
        attrs['type'] = type(self)
        attrs['bins'] = self.bins
        attrs['index_of_refraction'] = self.index_of_refraction
        attrs['distributionType'] = self.distributionType

        if 'layerbounderies' in dir(self):
            attrs['layerbounderies'] = self.layerbounderies

        storer.attrs.atmPy_attrs = attrs

        if 'data_fit_normal' in dir(self):
            table_name = table_name + '/data_fit_normal'
            hdf.put(table_name, self.data_fit_normal)
            storer = hdf.get_storer(table_name)
            storer.attrs.atmPy_attrs = None

        return hdf

    def zoom_diameter(self, start=None, end=None):
        sd = self.copy()
        if start:
            startIdx = array_tools.find_closest(sd.bins, start)
        else:
            startIdx = 0
        if end:
            endIdx = array_tools.find_closest(sd.bins, end)
        else:
            endIdx = len(self.bincenters)
        sd.data = self.data.iloc[:, startIdx:endIdx]
        sd.bins = self.bins[startIdx:endIdx + 1]
        self._update()
        return sd

    def _convert2otherDistribution(self, distType, verbose=False):
        self._mode_analysis = None
        return moments.convert(self, distType, verbose = verbose)

    def _get_mass_concentration(self):
        """'Mass concentration ($\mu g/m^{3}$)'"""
        dist = self.convert2dVdDp()
        volume_conc = dist.data * dist.binwidth

        vlc_all = volume_conc.sum(axis = 1) # nm^3/cm^3

        # if not self.properties.particle_density:
        #     raise ValueError('Please set the physical_property_density variable in g/cm^3')

        self.parameters4reductions._check_mixing_ratio_param_exist()

        if type(self.parameters4reductions.particle_density._dict['value']).__name__ in ['TimeSeries', 'VerticalProfile']:
            density = self.parameters4reductions.particle_density._dict['value'].copy()
            density = density.data['density']
        else:
            density = self.parameters4reductions.particle_density
            # density = self.physical_property_density #1.8 # g/cm^3

        density *= 1e-21 # g/nm^3
        mass_conc = vlc_all * density # g/cm^3
        mass_conc *= 1e6 # g/m^3
        mass_conc *= 1e6 # mug/m^3

        # mass_conc = self._normalize2temp_and_pressure(mass_conc)
        return mass_conc

    # def _normalize2temp_and_pressure(self, data):
    #     if
    #     data = atmosphere.normalize2pressure_and_temperature(data,P_is, P_shall, T_is, T_shall)
    #     return data

    def _get_mass_mixing_ratio(self):
        if not _panda_tools.ensure_column_exists(self.housekeeping.data, 'air_density_$g/m^3$', raise_error=False):
            _panda_tools.ensure_column_exists(self.housekeeping.data, 'temperature_K')
            _panda_tools.ensure_column_exists(self.housekeeping.data, 'pressure_Pa')
            gas = _gas_physics.Ideal_Gas_Classic()
            gas.temp = self.housekeeping.data['temperature_K']
            gas.pressure = self.housekeeping.data['pressure_Pa']
            self.housekeeping.data['air_density_$g/m^3$'] = gas.density_mass

        mass_conc = self.particle_mass_concentration.data.copy() #ug/m^3
        mass_conc = mass_conc.iloc[:,0]
        mass_conc *= 1e-6 #g/m^3
        mass_mix = mass_conc / (self.housekeeping.data['air_density_$g/m^3$'] - mass_conc)
        return mass_mix



    def _get_number_mixing_ratio(self):
        if not _panda_tools.ensure_column_exists(self.housekeeping.data, 'air_density_$number/m^3$', raise_error=False):
            _panda_tools.ensure_column_exists(self.housekeeping.data, 'temperature_K')
            _panda_tools.ensure_column_exists(self.housekeeping.data, 'pressure_Pa')
            gas = _gas_physics.Ideal_Gas_Classic()
            gas.temp = self.housekeeping.data['temperature_K']
            gas.pressure = self.housekeeping.data['pressure_Pa']
            self.housekeeping.data['air_density_$number/m^3$'] = gas.density_number


        numb_conc = self.particle_number_concentration.data.copy() # number/cc
        numb_conc = numb_conc.iloc[:,0]
        numb_conc *= 1e6 #number/m^3
        numb_mix = numb_conc/ (self.housekeeping.data['air_density_$number/m^3$'] - numb_conc)

        return numb_mix

    def _get_particle_concentration(self):
        """ Returns the sum of particles per line in data

        Returns
        -------
        int: if data has only one line
        pandas.DataFrame: else """
        sd = self.convert2numberconcentration()

        particles = sd.data.sum(axis = 1)

        # The code below is old and lead to problems when df contained NaNs
        # particles = _np.zeros(sd.data.shape[0])
        # for e, line in enumerate(sd.data.values):
        #     particles[e] = line.sum()
        if sd.data.shape[0] == 1:
            return particles[0]
        else:
            df = pd.DataFrame(particles,
                              # index=sd.data.index,
                              columns=['Particle number concentration #/$cm^3$'])

            # not sure what the purpose of the argument below was, however, it can
            # result in errors. If desired in future us: df = df.reindex(df.index.drop_duplicates())
            # df = df.drop_duplicates()
            return df

    def _get_surface_concentration(self):
        """ volume of particles per volume air"""

        sd = self.convert2dSdDp()

        surface_conc = sd.data * sd.binwidth

        sfc_all = surface_conc.sum(axis = 1) # nm^2/cm^3
        sfc_all = sfc_all * 1e-6 # um^2/cm^3
        label = 'Surface concentration $\mu m^2 / cm^{-3}$'
        sfc_df = pd.DataFrame(sfc_all, columns = [label])
        if type(self).__name__ == 'SizeDist':
            return sfc_df
        elif type(self).__name__ == 'SizeDist_TS':
            out =  _timeseries.TimeSeries(sfc_df)
            out._data_period = self._data_period
        elif type(self).__name__ == 'SizeDist_LS':
            out =  _vertical_profile.VerticalProfile(sfc_df)
        else:
            raise ValueError("Can't be! %s is not an option here"%(type(self).__name__))
        out._y_label = label
        return out



    def _get_volume_concentration(self):
        """ volume of particles per volume air"""

        sd = self.convert2dVdDp()

        volume_conc = sd.data * sd.binwidth

        vlc_all = volume_conc.sum(axis = 1) # nm^3/cm^3
        vlc_all = vlc_all * 1e-9 # um^3/cm^3
        vlc_df = pd.DataFrame(vlc_all, columns = ['volume concentration $\mu m^3 / cm^{3}$'])
        if type(self).__name__ == 'SizeDist':
            return  vlc_df
        elif type(self).__name__ == 'SizeDist_TS':
            out =  _timeseries.TimeSeries(vlc_df)
            out._data_period = self._data_period
        elif type(self).__name__ == 'SizeDist_LS':
            out = _vertical_profile.VerticalProfile(vlc_df)
        else:
            raise ValueError("Can't be! %s is not an option here"%(type(self).__name__))
        out._y_label = 'volume concentration $\mu m^3 / cm^{-3}$'
        return out


    def _hygro_growht_shift_data(self, data, bins, gf, ignore_data_nan = False):
        """data: 1D array
        bins: 1D array
        gf: float"""

        if _np.any(gf < 1):
            txt = 'Growth facotor smaller than 1 (is %s). Value adjusted to 1!!'%gf
            gf = 1.
            _warnings.warn(txt)

        # no shift if gf is 1
        if gf == 1.:
            out = {}
            out['bins'] = bins
            out['data'] = data
            out['num_extr_bins'] = 0
            return out
        # all data nan if gf is nan
        elif _np.isnan(gf):
            out = {}
            out['bins'] = bins
            data = data.copy()
            data[:] = _np.nan
            out['data'] = data
            out['num_extr_bins'] = 0
            return out
        # return as is if all data is nan
        elif _np.all(_np.isnan(data)):
            if not ignore_data_nan:
                out = {}
                out['bins'] = bins
                out['data'] = data
                out['num_extr_bins'] = 0
                return out


        bins_new = bins.copy()
        shifted = bins_new * gf

        ml = array_tools.find_closest(bins_new, shifted, how='closest_low')
        mh = array_tools.find_closest(bins_new, shifted, how='closest_high')

        if _np.any((mh - ml) > 1):
            raise ValueError('shifted bins spans over more than two of the original bins, programming required ;-)')

        no_extra_bins = bins_new[ml].shape[0] - bins_new[ml][bins_new[ml] != bins_new[ml][-1]].shape[0]

        ######### Ad bins to shift data into
        last_two = _np.log10(bins_new[- (no_extra_bins + 1):])
        step_width = last_two[-1] - last_two[-2]
        new_bins = _np.zeros(no_extra_bins)
        for i in range(no_extra_bins):
            new_bins[i] = _np.log10(bins_new[-1]) + ((i + 1) * step_width)
        newbins = 10**new_bins
        bins_new = _np.append(bins_new, newbins)
        shifted = (bins_new * gf)[:-no_extra_bins]

        ######## and again ########################

        ml = array_tools.find_closest(bins_new, shifted, how='closest_low')
        mh = array_tools.find_closest(bins_new, shifted, how='closest_high')

        if _np.any((mh - ml) > 1):
            raise ValueError('shifted bins spans over more than two of the original bins, programming required ;-)')


        ##### percentage of particles moved to next bin ...')

        shifted_w = shifted[1:] - shifted[:-1]

        fract_first = (bins_new[mh] - shifted)[:-1] / shifted_w
        fract_last = (shifted - bins_new[ml])[1:] / shifted_w

        data_new = _np.zeros(data.shape[0]+ no_extra_bins)
        data_new[no_extra_bins - 1:-1] += fract_first * data
        data_new[no_extra_bins:] += fract_last * data
        out = {}
        out['bins'] = bins_new
        out['data'] = data_new
        out['num_extr_bins'] = no_extra_bins
        return out



    def _update(self):
        """
        Resets properties so they are recalculated. This is usually necessary once you perform an operation on data.
        """
        self._optical_properties = None
        self._hygroscopicity = None
        self._mode_analysis = None
        self._uptodate_particle_number_concentration = False
        self._uptodate_particle_mass_concentration = False
        self._uptodate_particle_surface_concentration = False
        self._uptodate_particle_volume_concentration = False
        self._particle_mean_diameter = None
        self._submicron_volume_ratio = None
        self._normal_distribution_fits = None



class SizeDist_TS(SizeDist):
    """Returns a SizeDistribution_TS instance.

    Parameters:
    -----------
    data: pandas dataFrame with
         - column names (each name is something like this: '150-200')
         - index is time (at some point this should be arbitrary, convertable to altitude for example?)
    unit conventions:
         - diameters: nanometers
         - flowrates: cc (otherwise, axis label need to be adjusted an caution needs to be taken when dealing is AOD)
    distributionType:
         log normal: 'dNdlogDp','dSdlogDp','dVdlogDp'
         natural: 'dNdDp','dSdDp','dVdDp'
         number: 'dNdlogDp', 'dNdDp', 'numberConcentration'
         surface: 'dSdlogDp','dSdDp'
         volume: 'dVdlogDp','dVdDp'
       """
    def __init__(self,  *args, fill_data_gaps_with = None, ignore_data_gap_error = False, **kwargs):
        super(SizeDist_TS,self).__init__(*args,**kwargs)

        self._data_period = None
        # ignore_data_gap_error = False

        if type(fill_data_gaps_with).__name__ != 'NoneType':
            self.fill_gaps_with(what = fill_data_gaps_with)

        elif not ignore_data_gap_error:
            noofgaps = self.detect_gaps()
            if noofgaps > 0:
                raise ValueError(("There are {} gaps in the data. This might mean your instrument malfunctioned or there"
                                  "where time intervals where no particle was recorded. Use the fill_gaps_with function to"
                                  "fill data gaps with appropriate values or choose set the ignore_data_gap_error to True".format(noofgaps)))


        self.__particle_number_concentration = None
        self.__particle_mass_concentration = None
        self.__particle_mass_mixing_ratio = None
        self.__particle_number_mixing_ratio = None
        self._update()
        if not self.data.index.name:
            self.data.index.name = 'Time'

    close_gaps = _timeseries.close_gaps

    # def _update(self):
    #     self._uptodate_particle_number_concentration = False
    #     self._uptodate_particle_mass_concentration = False
    #     self._uptodate_particle_mass_mixing_ratio = False
    #     self._uptodate_particle_number_mixing_ratio = False
    #     self._uptodate_particle_surface_concentration = False
    #     self._uptodate_particle_volume_concentration = False
    #     self._optical_properties = None

    def detect_gaps(self, toleranz=1.95, return_all=False):
        idx = self.data.index
        dt = (idx[1:] - idx[:-1]) / _np.timedelta64(1, 's')

        med = _np.median(dt)
        mad = _robust.mad(dt)

        if mad == 0:
            noofgaps = 0
            dt = 0
            period = int(med)
        else:
            hist, edges = _np.histogram(dt[_np.logical_and((med - mad) < dt, dt < (med + mad))], bins=100)
            period = int(round((edges[hist.argmax()] + edges[hist.argmax() + 1]) / 2))
            noofgaps = dt[dt > toleranz * period].shape[0]

        if return_all:
                return {'index':idx,
                        'number of gaps':noofgaps,
                        'dt': dt,
                        'period (s)': period}
        else:
            return noofgaps

    def fill_gaps_with(self, what=0, toleranz=1.95):
        gaps = self.detect_gaps(toleranz=toleranz, return_all=True)
        idx = gaps['index']
        noofgaps = gaps['number of gaps']
        dt = gaps['dt']
        period = gaps['period (s)']
        print(type(toleranz),toleranz, type(period), period)
        for idxf, idxn, dtt in zip(idx[:-1][dt > toleranz * period], idx[1:][dt > toleranz * period],
                                   dt[dt > toleranz * period]):
            #     print(idxf, idxn, dtt)
            no2add = int(round(((idxn - idxf) / _np.timedelta64(1, 's')) / period)) - 1
            for i in range(no2add):
                newidx = idxf + _np.timedelta64((i + 1) * period, 's')
                self.data.loc[newidx, :] = what
        self.data.sort_index(inplace=True)
        return

    # def correct4ambient_LFE_tmp_difference(self):
    #     """corrects for temperature differences between ambient and instrument.
    #     The Problem is that the instrument samples at a constant flow at the temperature of
    #     the laminar flow element not the ambient temperatrue ... this corrects for it
    #     Make sure your housekeeping has a column named Temperature and one named Temperature_instrument"""
    #     sizedist = self.copy()
    #     if type(sizedist._correct4ambient_LFE_tmp_difference) != type(None):
    #         tcorr = (sizedist.housekeeping.data.Temperature_instrument + 273.15) / (sizedist.housekeeping.data.Temperature + 273.15)
    #         sizedist.data = sizedist.data.apply(lambda x: x * tcorr, axis=0)
    #         sizedist._update()
    #         sizedist._correct4ambient_LFE_tmp_difference = tcorr
    #     else:
    #         txt = "You can't apply this correction twice!"
    #         # raise ValueError(txt)
    #         _warnings.warn(txt)
    #     return sizedist

    # todo: declared deprecated on 2016-04-29
    def dprecated_convert2layerseries(self, hk, layer_thickness=10, force=False):
        """convertes the time series to a layer series.

        Note
        ----
        nan values are excluded when an average is taken over a the time that corresponds to the particular layer
        (altitude). If there are only nan values nan is returned and there is a gap in the Layerseries.

        The the housekeeping instance has to have a column called "Altitude" and which is monotonicly in- or decreasing

        Arguments
        ---------
        hk: housekeeping instance
        layer_thickness (optional): [10] thickness of each generated layer in meter"""
        _warnings.warn("convert2layerseries is deprecated and will be deleted in the future. Please use convert2verticalprofile (2016-04-29)")
        if any(_np.isnan(hk.data.Altitude)):
            txt = """The Altitude contains nan values. Either fix this first, eg. with pandas interpolate function"""
            raise ValueError(txt)

        if ((hk.data.Altitude.values[1:] - hk.data.Altitude.values[:-1]).min() < 0) and (
                    (hk.data.Altitude.values[1:] - hk.data.Altitude.values[:-1]).max() > 0):
            if force:
                hk.data = hk.data.sort(columns='Altitude')
            else:
                txt = ('Given altitude data is not monotonic. This is not possible (yet). Use force if you\n'
                       'know what you are doing')
                raise ValueError(txt)

        start_h = round(hk.data.Altitude.values.min() / layer_thickness) * layer_thickness #what did I do here?!?!
        end_h = round(hk.data.Altitude.values.max() / layer_thickness) * layer_thickness

        layer_edges = _np.arange(start_h, end_h, layer_thickness)
        empty_frame = pd.DataFrame(columns=self.data.columns)
        lays = SizeDist_LS(empty_frame, self.bins, self.distributionType, None)

        for e, end_h_l in enumerate(layer_edges[1:]):
            start_h_l = layer_edges[e]
            layer = hk.data.Altitude.iloc[
                _np.where(_np.logical_and(start_h_l < hk.data.Altitude.values, hk.data.Altitude.values < end_h_l))]
            start_t = layer.index.min()
            end_t = layer.index.max()
            dist_tmp = self.zoom_time(start=start_t, end=end_t)
            avrg = dist_tmp.average_overAllTime()
            # return avrg,lays
            lays.add_layer(avrg, (start_h_l, end_h_l))
        lays.parent_dist_TS = self
        lays.parent_timeseries = hk

        data = hk.data.copy()
        data['Time_UTC'] = data.index
        data.index = data.Altitude
        data = data.sort_index()
        if not data.index.is_unique:  # this is needed in case there are duplicate indeces
            grouped = data.groupby(level=0)
            data = grouped.last()

        # lays.housekeeping = data
        data = data.reindex(lays.layercenters, method='nearest')
        lays.housekeeping = _vertical_profile.VerticalProfile(data)
        return lays

    def convert2verticalprofile(self, layer_thickness=2):
        sd = self.copy()
        # if not np.array_equal(sd.data.index.values, sd.housekeeping.data.index.values): raise ValueError()
        # if not np.array_equal(sd.data.index.values, sd.housekeeping.data.index.values): raise ValueError()

        sd.data.index = sd.housekeeping.data.Altitude
        sd.data.sort_index(inplace=True)

        sd.housekeeping.data.index = sd.housekeeping.data.Altitude
        sd.housekeeping.data.sort_index(inplace=True)

        start = _np.floor(sd.housekeeping.data.Altitude.min())
        end = _np.ceil(sd.housekeeping.data.Altitude.max())

        layerbounderies = _np.arange(start, end + 1, layer_thickness)
        layerbounderies = _np.array([layerbounderies[0:-1], layerbounderies[1:]]).transpose()

        index = _np.apply_along_axis(lambda x: x.sum(), 1, layerbounderies) / 2.
        df = pd.DataFrame(_np.zeros((layerbounderies.shape[0], sd.data.shape[1])), index=index, columns=sd.data.columns)

        dfhk = pd.DataFrame(_np.zeros((layerbounderies.shape[0], sd.housekeeping.data.shape[1])), index=index, columns=sd.housekeeping.data.columns)

        for l in layerbounderies:
            where = _np.where(_np.logical_and(sd.data.index > l[0], sd.data.index < l[1]))[0]
            mean = sd.data.iloc[where, :].mean()
            df.loc[(l[0] + l[1]) / 2] = mean

            mean = sd.housekeeping.data.iloc[where, :].mean()
            dfhk.loc[(l[0] + l[1]) / 2] = mean

        dist_ls = SizeDist_LS(df, sd.bins, sd.distributionType, layerbounderies)
        dist_ls.housekeeping = _vertical_profile.VerticalProfile(dfhk)
        return dist_ls

    def fit_normal_deprecated(self, log=True, p0=[10, 180, 0.2]):
        """ Fits a single normal distribution to each line in the data frame.

        Returns
        -------
        pandas DataFrame instance (also added to namespace as data_fit_normal)

        """

        super(SizeDist_TS, self).fit_normal(log=log, p0=p0)
        self.data_fit_normal.index = self.data.index
        return self.data_fit_normal


    def deprecated_apply_hygro_growth(self, kappa, RH = None, how='shift_data', adjust_refractive_index = True):
        """ see docstring of atmPy.sizedistribution.SizeDist for more information
        Parameters
        ----------
        kappa: float
        RH: bool, float, or array.
            If None, RH from self.housekeeping will be taken"""
        # import pdb; pdb.set_trace()
        if not _np.any(RH):
            pandas_tools.ensure_column_exists(self.housekeeping.data, 'Relative_humidity')
            RH = self.housekeeping.data.Relative_humidity.values
        # return kappa,RH,how
        sd = super(SizeDist_TS,self).deprecated_apply_hygro_growth(kappa, RH, how = how, adjust_refractive_index = adjust_refractive_index)
        # sd = out['size_distribution']
        # gf = out['growth_factor']
        sd_TS = SizeDist_TS(sd.data, sd.bins, sd.distributionType, fixGaps=False)
        sd_TS.parameters4reductions.refractive_index = sd.parameters4reductions.refractive_index
        # sd_TS._SizeDist_growth_factor = sd.growth_factor
        # out['size_distribution'] = sd_LS
        return sd_TS

    def deprecated_apply_growth(self, growth_factor, how='shift_data'):
        """ see docstring of atmPy.sizedistribution.SizeDist for more information
        Parameters
        ----------
        kappa: float
        RH: bool, float, or array.
            If None, RH from self.housekeeping will be taken"""

        sd = super(SizeDist_TS,self).deprecated_apply_growth(growth_factor, how = how)
        # sd = out['size_distribution']
        # gf = out['growth_factor']
        sd_TS = SizeDist_TS(sd.data, sd.bins, sd.distributionType, fixGaps=False)
        sd_TS.parameters4reductions.refractive_index = sd.parameters4reductions.refractive_index
        # sd_TS._SizeDist_growth_factor = sd.growth_factor
        sd_TS._data_period = self._data_period
        # out['size_distribution'] = sd_LS
        return sd_TS

    # def calculate_optical_properties(self, wavelength, n = None, noOfAngles=100):
    #     # opt = super(SizeDist_TS,self).calculate_optical_properties(wavelength, n = None, AOD = False, noOfAngles=100)
    #     if not _np.any(n):
    #         n = self.index_of_refraction
    #     if not _np.any(n):
    #         txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #         raise ValueError(txt)
    #
    #     out = optical_properties.size_dist2optical_properties(self, wavelength, n,
    #                                                           aod=False,
    #                                                           noOfAngles=noOfAngles)
    #     # opt_properties = optical_properties.OpticalProperties(out, self.bins)
    #     # opt._data_period = self._data_period
    #     return out

    def _getXYZ(self):
        """
        This will create three arrays, so when plotted with pcolor each pixel will represent the exact bin width
        """
        binArray = _np.repeat(_np.array([self.bins]), self.data.index.shape[0], axis=0)
        timeArray = _np.repeat(_np.array([self.data.index.values]), self.bins.shape[0], axis=0).transpose()
        ext = _np.array([_np.zeros(self.data.index.values.shape)]).transpose()
        Z = _np.append(self.data.values, ext, axis=1)
        return timeArray, binArray, Z

    def get_timespan(self):
        return self.data.index.min(), self.data.index.max()

    # TODO: Fix plot options such as showMinorTickLabels
    def plot(self,
             vmax=None,
             vmin=None,
             norm='linear',
             showMinorTickLabels=True,
             # removeTickLabels=["700", "900"],
             ax=None,
             fit_pos=False,
             cmap=plt_tools.get_colorMap_intensity(),
             colorbar=True):

        """ plots an intensity plot of all data

        Arguments
        ---------
        scale (optional): ('log',['linear']) - defines how the z-direction is scaled
        vmax
        vmin
        show_minor_tickLabels:
        cma:
        fit_pos: bool[True]. Optional
            plots the position of a fitted normal distribution onto the plot.
            in order for this to work execute fit_normal
        ax (optional):  axes instance [None] - option to plot on existing axes

        Returns
        -------
        f,a,pc,cb (figure, axis, pcolormeshInstance, colorbar)

        """
        X, Y, Z = self._getXYZ()
        Z = _np.ma.masked_invalid(Z)

        if type(ax).__name__ in _axes_types:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()
            f.autofmt_xdate()

        if norm == 'log':
            norm = LogNorm()
        elif norm == 'linear':
            norm = None

        pc = a.pcolormesh(X, Y, Z, vmin=vmin, vmax=vmax, norm=norm, cmap=cmap)
        a.set_yscale('log')
        a.set_ylim((self.bins[0], self.bins[-1]))
        a.set_xlabel('Time (UTC)')

        a.get_yaxis().set_tick_params(direction='out', which='both')
        a.get_xaxis().set_tick_params(direction='out', which='both')

        if self.distributionType == 'calibration':
            a.set_ylabel('Amplitude (digitizer bins)')
        else:
            a.set_ylabel('Diameter (nm)')
        if colorbar:
            cb = f.colorbar(pc)
            label = get_label(self.distributionType)
            cb.set_label(label)
        else:
            cb = get_label(self.distributionType)

        # if self.distributionType != 'calibration':
            # a.yaxis.set_major_formatter(plt.FormatStrFormatter("%i"))

            # f.canvas.draw()  # this is important, otherwise the ticks (at least in case of minor ticks) are not created yet
        if showMinorTickLabels:
            minf = plt_tools.get_formatter_minor_log()
            a.yaxis.set_minor_formatter(minf)
                # a.yaxis.set_minor_formatter(plt.FormatStrFormatter("%i"))
                # ticks = a.yaxis.get_minor_ticks()
                # for i in ticks:
                #     if i.label.get_text() in removeTickLabels:
                #         i.label.set_visible(False)

        if fit_pos:
            # if 'data_fit_normal' in dir(self):
            f,a,_ = self.normal_distribution_fits.plot_fitres(show_amplitude=False, show_width=False, ax = a)
            g = a.get_lines()[-1]
            g.set_color('m')
            g.set_label('normal dist. center')
            # a.plot(self.data.index, self.data_fit_normal.Pos, color='m', linewidth=2, label='normal dist. center')
            leg = a.legend(loc = 1, fancybox=True, framealpha=0.5)
            leg.draw_frame(True)

        return f, a, pc, cb

    # def plot_fitres(self):
    #     """ Plots the results from fit_normal"""
    #
    #     f, a = plt.subplots()
    #     data = self.data_fit_normal.dropna()
    #     a.fill_between(data.index, data.Sigma_high, data.Sigma_low,
    #                    color=plt_tools.color_cycle[0],
    #                    alpha=0.5,
    #                    )
    #     a.plot(data.index.values, data.Pos.values, color=plt_tools.color_cycle[0], linewidth=2, label='center')
    #     # data.Pos.plot(ax=a, color=plt_tools.color_cycle[0], linewidth=2, label='center')
    #     a.legend(loc=2)
    #     a.set_ylabel('Particle diameter (nm)')
    #     a.set_xlabel('Altitude (m)')
    #
    #     a2 = a.twinx()
    #     # data.Amp.plot(ax=a2, color=plt_tools.color_cycle[1], linewidth=2, label='amplitude')
    #     a2.plot(data.index.values, data.Amp.values, color=plt_tools.color_cycle[1], linewidth=2, label='amplitude')
    #     a2.legend()
    #     a2.set_ylabel('Amplitude - %s' % (get_label(self.distributionType)))
    #     f.autofmt_xdate()
    #     return f, a, a2

    # def plot_particle_concentration(self, ax=None, label=None):
    #     """Plots the particle rate as a function of time.
    #
    #     Parameters
    #     ----------
    #     ax: matplotlib.axes instance, optional
    #         perform plot on these axes.
    #
    #     Returns
    #     -------
    #     matplotlib.axes instance
    #
    #     """
    #
    #     if type(ax).__name__ in _axes_types:
    #         color = plt_tools.color_cycle[len(ax.get_lines())]
    #         f = ax.get_figure()
    #     else:
    #         f, ax = plt.subplots()
    #         color = plt_tools.color_cycle[0]
    #
    #     # layers = self.convert2numberconcentration()
    #
    #     particles = self.get_particle_concentration().dropna()
    #
    #     ax.plot(particles.index.values, particles.Count_rate.values, color=color, linewidth=2)
    #
    #     if label:
    #         ax.get_lines()[-1].set_label(label)
    #         ax.legend()
    #
    #     ax.set_xlabel('Time (UTC)')
    #     ax.set_ylabel('Particle number concentration (cm$^{-3})$')
    #     if particles.index.dtype.type.__name__ == 'datetime64':
    #         f.autofmt_xdate()
    #     return ax

    def zoom_time(self, start=None, end=None):
        """
        2014-11-24 16:02:30
        """
        dist = self.copy()
        dist.data = dist.data.truncate(before=start, after=end)
        if dist.housekeeping:
            dist.housekeeping = self.housekeeping.zoom_time(start=start, end=end)

        dist._update()
        return dist


    def average_time(self, window='10s'):
        """returns a copy of the sizedistribution_TS with reduced size by 
        averaging over a given window. Thesedays it is merely a wrap around
        pandas resample function

        Arguments
        ---------
        window: str
            see pandas resample doc

        Returns
        -------
        SizeDistribution_TS instance
            copy of current instance with resampled data frame
        """

        dist = self.copy()
        dist.data = dist.data.resample(window, label='left').mean()
        # dist._data_period = _np.timedelta64(window[0], window[1]) / _np.timedelta64(1, 's')
        dist._data_period = int(pd.to_timedelta(window)/pd.to_timedelta('1s'))
        if dist.distributionType == 'calibration':
            dist.data.values[_np.where(_np.isnan(self.data.values))] = 0

        if dist.housekeeping:
            dist.housekeeping = self.housekeeping.average_time(window = window)


        dist._update()
        return dist

    def average_overAllTime(self, sigma = 0, minmax = False, percentile = False):
        """
        averages over the entire dataFrame and returns a single sizedistribution (numpy.ndarray)

        Args
        ----
        minmax: bool
            returns in addition size distribution that represent the min and max values
        percentile: float
            percentile to be calculated and returned in addition to averaged sizedist
        sigma (!experimental!): int
            if not ==0 this function will additionally return the std according to sigma
        """
        singleHist = _np.zeros(self.data.shape[1])
        if sigma:
            singleHist_std = _np.zeros(self.data.shape[1])

        elif minmax:
            singleHist_min = _np.zeros(self.data.shape[1])
            singleHist_max = _np.zeros(self.data.shape[1])

        elif percentile:
            # if type(percentile).__name__ in (int, float):
            #     percentile = [percentile]
            singleHist_percentile = _np.zeros(self.data.shape[1])

        for i in range(self.data.shape[1]):
            line = self.data.values[:, i]
            singleHist[i] = _np.average(line[~_np.isnan(line)])
            if sigma:
                singleHist_std[i] = _np.std(line[~_np.isnan(line)])
            elif minmax:
                try:
                    singleHist_min[i] = _np.min(line[~_np.isnan(line)])
                    singleHist_max[i] = _np.max(line[~_np.isnan(line)])
                except ValueError:
                    pass
            elif percentile:
                try:
                    singleHist_percentile[i] = _np.percentile(line[~_np.isnan(line)], percentile)
                except IndexError:
                    pass

        data = pd.DataFrame(_np.array([singleHist]), columns=self.data.columns)
        avgDist = SizeDist(data, self.bins, self.distributionType)

        if sigma:
            data_low = pd.DataFrame(_np.array([singleHist - singleHist_std]), columns=self.data.columns)
            std_dist_low = SizeDist(data_low, self.bins, self.distributionType)
            data_high = pd.DataFrame(_np.array([singleHist + singleHist_std]), columns=self.data.columns)
            std_dist_high = SizeDist(data_high, self.bins, self.distributionType)
            return avgDist, std_dist_low, std_dist_high
        elif minmax:
            data = pd.DataFrame(_np.array([singleHist_min]), columns=self.data.columns)
            mindist = SizeDist(data, self.bins, self.distributionType)
            data = pd.DataFrame(_np.array([singleHist_max]), columns=self.data.columns)
            maxdist = SizeDist(data, self.bins, self.distributionType)
            return avgDist, mindist, maxdist

        elif percentile:
            data = pd.DataFrame(_np.array([singleHist_percentile]), columns=self.data.columns)
            percDist = SizeDist(data, self.bins, self.distributionType)
            return avgDist, percDist

        else:
            return avgDist

    @property
    def particle_number_concentration(self):
        if not _np.any(self.__particle_number_concentration) or not self._uptodate_particle_number_concentration:
            self.__particle_number_concentration = _timeseries.TimeSeries(self._get_particle_concentration())
            self.__particle_number_concentration._y_label = 'Particle number concentration #/$cm^3$'
            self.__particle_number_concentration._x_label = 'Time'
            self._uptodate_particle_number_concentration = True
            self.__particle_number_concentration._data_period = self._data_period
        return self.__particle_number_concentration

    @property
    def particle_mean_diameter(self):
        sup = super().particle_mean_diameter
        return _timeseries.TimeSeries(sup, sampling_period = self._data_period)

    @property
    def particle_mass_concentration(self):
        if not _np.any(self.__particle_mass_concentration) or not self._uptodate_particle_mass_concentration:
            mass_conc = self._get_mass_concentration()
            mass_conc = pd.DataFrame(mass_conc, columns = ['Mass concentration ($\mu g/m^{3}$)'])
            self.__particle_mass_concentration = _timeseries.TimeSeries(mass_conc)
            self.__particle_mass_concentration._y_label = 'Mass concentration ($\mu g/m^{3}$)'
            self.__particle_mass_concentration._x_label =  'Time'
            self.__particle_mass_concentration._data_period = self._data_period
            self._uptodate_particle_mass_concentration = True
        return self.__particle_mass_concentration


    @property
    def particle_mass_mixing_ratio(self):
        if not _np.any(self.__particle_mass_mixing_ratio) or not self._uptodate_particle_mass_mixing_ratio:
            mass_mix = self._get_mass_mixing_ratio()
            ylabel = 'Particle mass mixing ratio'
            mass_mix = pd.DataFrame(mass_mix)
            self.__particle_mass_mixing_ratio = _timeseries.TimeSeries(mass_mix)
            self.__particle_mass_mixing_ratio._data_period = self._data_period
            self.__particle_mass_mixing_ratio._y_label = ylabel
            self.__particle_mass_mixing_ratio._x_label = 'Time'
            self._uptodate_particle_mass_mixing_ratio = True
        return self.__particle_mass_mixing_ratio

    @property
    def particle_number_mixing_ratio(self):
        if not _np.any(self.__particle_number_mixing_ratio) or not self._uptodate_particle_number_mixing_ratio:
            number_mix = self._get_number_mixing_ratio()
            ylabel = 'Particle number mixing ratio'
            number_mix = pd.DataFrame(number_mix)
            self.__particle_number_mixing_ratio = _timeseries.TimeSeries(number_mix)
            self.__particle_number_mixing_ratio._data_period = self._data_period
            self.__particle_number_mixing_ratio._y_label = ylabel
            self.__particle_number_mixing_ratio._x_label = 'Time'
            self._uptodate_particle_number_mixing_ratio = True
        return self.__particle_number_mixing_ratio

    @property
    def optical_properties(self):
        if not self._optical_properties:
            self._optical_properties = optical_properties.OpticalProperties_TS(self)
        return self._optical_properties

    # @property
    # def optical_properties(self):
    #     self.__optical_properties = optical_properties.OpticalProperties_TS(self)
    #     return self.__optical_properties
    #
    # @optical_properties.setter
    # def optical_properties(self, value):
    #     self.__optical_properties = value


class SizeDist_LS(SizeDist):
    """
    Parameters
    ----------
    data:    pandas DataFrame ...
    bins:    array
    distributionType:   str
    layerbounderies:    array shape(n_layers,2)

    OLD
    ---
    data: pandas dataFrame with
                 - column names (each name is something like this: '150-200')
                 - altitude (at some point this should be arbitrary, convertable to altitude for example?)
       unit conventions:
             - diameters: nanometers
             - flowrates: cc (otherwise, axis label need to be adjusted an caution needs to be taken when dealing is AOD) 
       distributionType:  
             log normal: 'dNdlogDp','dSdlogDp','dVdlogDp'
             natural: 'dNdDp','dSdDp','dVdDp'
             number: 'dNdlogDp', 'dNdDp', 'numberConcentration'
             surface: 'dSdlogDp','dSdDp'
             volume: 'dVdlogDp','dVdDp'

    """

    def __init__(self, data, bins, distributionType, layerbounderies):
        super(SizeDist_LS, self).__init__(data, bins, distributionType,
                                          # fixGaps=False
                                          )
        if type(layerbounderies).__name__ == 'NoneType':
            self.layerbounderies = _np.empty((0, 2))
            # self.layercenters = _np.array([])
        else:
            self.layerbounderies = layerbounderies

        self.__particle_number_concentration = None
        self.__particle_mass_concentration = None
        self.__particle_mass_mixing_ratio = None
        self.__particle_number_mixing_ratio = None
        self._optical_properties = None
        # self.sup_optical_properties_wavelength = None
        self._SizeDist__sup_opt_aod = True
        self._update()

    # todo: I have the impression this is redundant. I think the reason I introduced this, was that I did not understand how inheritance of properties work .. uuuuh what a mess :-)
    # def _update(self):
    #     self._uptodate_particle_number_concentration = False
    #     self._uptodate_particle_mass_concentration = False
    #     self._uptodate_particle_mass_mixing_ratio = False
    #     self._uptodate_particle_number_mixing_ratio = False

    @property
    def optical_properties(self):
        if not self._optical_properties:
            self._optical_properties = optical_properties.OpticalProperties_VP(self)
        return self._optical_properties

    # @property
    # def optical_properties(self):
    #     self.__optical_properties = optical_properties.OpticalProperties_VP(self)
    #     return self.__optical_properties
    #
    # @optical_properties.setter
    # def optical_properties(self, value):
    #     self.__optical_properties = value

    @property
    def housekeeping(self):
        return self._housekeeping

    @housekeeping.setter
    def housekeeping(self, value):
        self._housekeeping = value #had to overwrite the setter since there is no align programmed yet for verticle profiles ... do it!

    @property
    def layercenters(self):
        return self.__layercenters

    @property
    def layerbounderies(self):
        return self.__layerbouderies

    @layerbounderies.setter
    def layerbounderies(self,lb):
        self.__layerbouderies = lb
        # newlb = _np.unique(self.layerbounderies.flatten()) # the unique is sorting the data, which is not reallyt what we want!
        # self.__layercenters = (newlb[1:] + newlb[:-1]) / 2.
        self.__layercenters = (self.layerbounderies[:,0] + self.layerbounderies[:,1]) / 2.
        self.data.index = self.layercenters


    @property
    def particle_number_concentration(self):
        if not _np.any(self.__particle_number_concentration) or not self._uptodate_particle_number_concentration:
            self.__particle_number_concentration = _vertical_profile.VerticalProfile(self._get_particle_concentration())
            self.__particle_number_concentration._x_label = 'Particle number concentration (#/$cm^3#)'
            self._uptodate_particle_number_concentration = True
        return self.__particle_number_concentration

    @property
    def particle_mass_concentration(self):
        if not _np.any(self.__particle_mass_concentration) or not self._uptodate_particle_mass_concentration:
            mass_conc = self._get_mass_concentration()
            mass_conc = pd.DataFrame(mass_conc, columns = ['Mass concentration ($\mu g/m^{3}$)'])
            self.__particle_mass_concentration = _vertical_profile.VerticalProfile(mass_conc)
            self.__particle_mass_concentration._x_label = 'Mass concentration ($\mu g/m^{3}$)'
            self.__particle_mass_concentration._y_label =  'Altitde'

            self._uptodate_particle_mass_concentration = True
        return self.__particle_mass_concentration

    @property
    def particle_mass_mixing_ratio(self):
        if not _np.any(self.__particle_mass_mixing_ratio) or not self._uptodate_particle_mass_mixing_ratio:
            mass_mix = self._get_mass_mixing_ratio()
            ylabel = 'Particle mass mixing ratio'
            mass_mix = pd.DataFrame(mass_mix, columns = [ylabel])
            self.__particle_mass_mixing_ratio = _vertical_profile.VerticalProfile(mass_mix)

            self.__particle_mass_mixing_ratio.data.index.name = 'Altitude'
            self.__particle_mass_mixing_ratio._x_label = ylabel
            self.__particle_mass_mixing_ratio._y_label = 'Altitude'
            self._uptodate_particle_mass_mixing_ratio = True
        return self.__particle_mass_mixing_ratio

    @property
    def particle_number_mixing_ratio(self):
        if not _np.any(self.__particle_number_mixing_ratio) or not self._uptodate_particle_number_mixing_ratio:
            number_mix = self._get_number_mixing_ratio()
            ylabel = 'Particle number mixing ratio'
            number_mix = pd.DataFrame(number_mix, columns = [ylabel])
            self.__particle_number_mixing_ratio = _vertical_profile.VerticalProfile(number_mix)
            self.__particle_number_mixing_ratio.data.index.name = 'Altitude'
            self.__particle_number_mixing_ratio._x_label = ylabel
            self.__particle_number_mixing_ratio._y_label = 'Altitude'
            self._uptodate_particle_number_mixing_ratio = True
        return self.__particle_number_mixing_ratio

    # @property
    # def optical_properties(self):
    #     if not self.__optical_properties:
    #         if not _np.any(self.index_of_refraction):
    #             txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #             raise ValueError(txt)
    #         if not self.sup_optical_properties_wavelength:
    #             txt = 'Please provied wavelength by setting the attribute sup_optical_properties_wavelength (in nm).'
    #             raise AttributeError(txt)
    #         self.__optical_properties = optical_properties.size_dist2optical_properties(self, self.sup_optical_properties_wavelength, self.index_of_refraction, aod = True, noOfAngles=100)
    #         # opt_properties = optical_properties.OpticalProperties(out, parent = self)
    #         # opt_properties.wavelength = wavelength #should be set in OpticalProperty class
    #         # opt_properties.index_of_refractio = n
    #         #opt_properties.angular_scatt_func = out['angular_scatt_func']  # This is the formaer phase_fct, but since it is the angular scattering intensity, i changed the name
    #         # opt_properties.parent_dist = self
    #     return self.__optical_properties

    def deprecated_apply_hygro_growth(self, kappa, RH = None, how='shift_data'):
        """ see docstring of atmPy.sizedistribution.SizeDist for more information
        Parameters
        ----------
        kappa: float
        RH: bool, float, or array.
            If None, RH from self.housekeeping will be taken"""
        if not _np.any(RH):
            pandas_tools.ensure_column_exists(self.housekeeping.data, 'Relative_humidity')
            RH = self.housekeeping.data.Relative_humidity.values
        # return kappa,RH,how
        sd = super(SizeDist_LS,self).deprecated_apply_hygro_growth(kappa, RH, how = how)
        # sd = out['size_distribution']
        # gf = out['growth_factor']
        sd_LS = SizeDist_LS(sd.data, sd.bins, sd.distributionType, self.layerbounderies,
                            # fixGaps=False
                            )
        # sd_LS.hygroscopic_growth_settings.refractive_index = sd.hygroscopic_growth_settings.refractive_index
        sd_LS.parameters4reductions.refractive_index = sd.parameters4reductions.refractive_index
        sd_LS.parameters4reductions.wavelength = sd.parameters4reductions.wavelength
        # sd_LS._SizeDist_growth_factor = sd.growth_factor
        # out['size_distribution'] = sd_LS
        return sd_LS

    def deprecated_apply_growth(self, growth_factor, how='shift_data'):
        """ see docstring of atmPy.sizedistribution.SizeDist for more information
        Parameters
        ----------
        kappa: float
        RH: bool, float, or array.
            If None, RH from self.housekeeping will be taken"""

        sd = super(SizeDist_LS,self).deprecated_apply_growth(growth_factor, how = how)
        # sd = out['size_distribution']
        # gf = out['growth_factor']
        sd_LS = SizeDist_LS(sd.data, sd.bins, sd.distributionType, self.layerbounderies,
                            # fixGaps=False
                            )

        sd_LS.parameters4reductions.refractive_index = sd.parameters4reductions.refractive_index
        # pdb.set_trace()
        sd_LS.parameters4reductions.wavelength = sd.parameters4reductions.wavelength
        # sd_LS._SizeDist_growth_factor = sd.growth_factor
        # out['size_distribution'] = sd_LS
        return sd_LS

    def deprecated_calculate_angstromex(self, wavelengths=[460.3, 550.4, 671.2, 860.7], n=1.455):
        """Calculates the Anstrome coefficience (overall, layerdependent)

        Parameters
        ----------
        wavelengths:    array-like, optional.
            the angstrom coefficient will be calculated based on the AOD of these wavelength values (in nm)
        n:              float, optional.
            index of refraction used in the underlying mie calculation.

        Returns
        -------
        Angstrom exponent, float
        List containing the OpticalProperties instances for the different wavelengths

        New Attributes
        --------------
        angstromexp:        float
            the resulting angstrom exponent
        angstromexp_fit:    pandas instance.
            AOD and fit result as a function of wavelength
        angstromexp_LS:     pandas instance.
            angstrom exponent as a function of altitude
        """

        AOD_list = []
        AOD_dict = {}
        for w in wavelengths:
            AOD = self.calculate_optical_properties(w, n)  # calculate_AOD(wavelength=w, n=n)
            #     opt= sizedistribution.OpticalProperties(AOD, dist_LS.bins)
            AOD_list.append({'wavelength': w, 'opt_inst': AOD})
            AOD_dict['%.1f' % w] = AOD

        eg = AOD_dict[list(AOD_dict.keys())[0]]

        wls = AOD_dict.keys()
        wls_a = _np.array(list(AOD_dict.keys())).astype(float)
        ang_exp = []
        ang_exp_std = []
        ang_exp_r_value = []
        for e, el in enumerate(eg.layercenters):
            AODs = _np.array([AOD_dict[wl].data_orig['AOD_layer'].values[e][0] for wl in wls])
            slope, intercept, r_value, p_value, std_err = stats.linregress(_np.log10(wls_a), _np.log10(AODs))
            ang_exp.append(-slope)
            ang_exp_std.append(std_err)
            ang_exp_r_value.append(r_value)
        # break
        ang_exp = _np.array(ang_exp)
        ang_exp_std = _np.array(ang_exp_std)
        ang_exp_r_value = _np.array(ang_exp_r_value)

        tmp = _np.array([[float(i), AOD_dict[i].AOD] for i in AOD_dict.keys()])
        wavelength, AOD = tmp[_np.argsort(tmp[:, 0])].transpose()
        slope, intercept, r_value, p_value, std_err = stats.linregress(_np.log10(wavelength), _np.log10(AOD))

        self.angstromexp = -slope
        aod_fit = _np.log10(wavelengths) * slope + intercept
        self.angstromexp_fit = pd.DataFrame(_np.array([AOD, 10 ** aod_fit]).transpose(), index=wavelength,
                                            columns=['data', 'fit'])

        self.angstromexp_LS = pd.DataFrame(_np.array([ang_exp, ang_exp_std, ang_exp_r_value]).transpose(),
                                           index=self.layercenters,
                                           columns=['ang_exp', 'standard_dif', 'correlation_coef'])
        self.angstromexp_LS.index.name = 'layercenter'

        return -slope, AOD_dict

    # # todo: this function appears multiple times, can easily be inherited ... probably not anymor, still testing
    # def calculate_optical_properties(self, wavelength, n = None, noOfAngles=100):
    #     if not n:
    #         n = self.index_of_refraction
    #     if not n:
    #         txt = 'Refractive index is not specified. Either set self.index_of_refraction or set optional parameter n.'
    #         raise ValueError(txt)
    #     out = _calculate_optical_properties(self, wavelength, n, aod = True, noOfAngles=noOfAngles)
    #     opt_properties = OpticalProperties(out, self.bins)
    #     opt_properties.wavelength = wavelength
    #     opt_properties.index_of_refractio = n
    #     opt_properties.angular_scatt_func = out['angular_scatt_func']  # This is the formaer phase_fct, but since it is the angular scattering intensity, i changed the name
    #     opt_properties.parent_dist_LS = self
    #     return opt_properties

    # def calculate_optical_properties(self, wavelength, n = None, AOD = True, noOfAngles=100):
    #     opt = super(SizeDist_LS,self).calculate_optical_properties(wavelength, n = None, AOD = True, noOfAngles=100)
    #     return opt


    def add_layer(self, sd, layerboundery):
        """
        Adds a sizedistribution instance to the layerseries.
        layerboundery

        Parameters
        ----------
        sd:
        layerboundary:
        """
        if len(layerboundery) != 2:
            raise ValueError('layerboundery has to be of length 2')

        sd = sd._convert2otherDistribution(self.distributionType)

        layerbounderies = _np.append(self.layerbounderies, _np.array([layerboundery]), axis=0)

        layerbounderiesU = _np.unique(layerbounderies)

        if (_np.where(layerbounderiesU == layerboundery[1])[0] - _np.where(layerbounderiesU == layerboundery[0])[0])[
            0] != 1:
            raise ValueError('The new layer is overlapping with an existing layer!')
        self.data = self.data.append(sd.data)
        self.layerbounderies = layerbounderies
        # self.layerbounderies.sort(axis=0)
        #
        # layercenter = _np.array(layerboundery).sum() / 2.
        # self.layercenters = _np.append(self.layercenters, layercenter)
        # self.layercenters.sort()
        # sd.data.index = _np.array([layercenter])

        # self.data = self.data.append(sd.data)
        return

    def _getXYZ(self):
        """
        This will create three arrays, so when plotted with pcolor each pixel will represent the exact bin width
        """
        binArray = _np.repeat(_np.array([self.bins]), self.data.index.shape[0], axis=0)
        layerArray = _np.repeat(_np.array([self.data.index.values]), self.bins.shape[0], axis=0).transpose()
        ext = _np.array([_np.zeros(self.data.index.values.shape)]).transpose()
        Z = _np.append(self.data.values, ext, axis=1)
        return layerArray, binArray, Z

    def plot_eachLayer(self, a=None, normalize=False):
        """
        Plots the distribution of each layer in one plot.

        Returns
        -------
        Handles to the figure and axes of the plot
        """
        if not a:
            f, a = _plt.subplots()
        else:
            f = None
            pass
        for iv in self.data.index.values:
            if normalize:
                a.plot(self.bincenters, self.data.loc[iv, :] / self.data.loc[iv, :].max(), label='%i' % iv)
            else:
                a.plot(self.bincenters, self.data.loc[iv, :], label='%i' % iv)
        a.set_xlabel('Particle diameter (nm)')
        a.set_ylabel(get_label(self.distributionType))
        a.legend()
        a.semilogx()
        return f, a


    def plot(self, vmax=None, vmin=None, scale='linear', show_minor_tickLabels=True,
             removeTickLabels=["500", "700", "800", "900"],
             plotOnTheseAxes=False,
             cmap=plt_tools.get_colorMap_intensity(),
             fit_pos=False,
             ax=None,
             colorbar = True):
        """ plots and returns f,a,pc,cb (figure, axis, pcolormeshInstance, colorbar)

        Arguments
        ---------
        scale (optional): ('log',['linear']) - defines how the z-direction is scaled
        vmax
        vmin
        show_minor_tickLabels:
        cma:
        fit_pos (optional): bool [True] - plots the position of a fitted normal distribution onto the plot.
                                          in order for this to work execute fit_normal
        ax (optional):  axes instance [None] - option to plot on existing axes

        """
        X, Y, Z = self._getXYZ()
        Z = _np.ma.masked_invalid(Z)
        if type(ax).__name__ in _axes_types:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()
            # f.autofmt_xdate()

        if scale == 'log':
            scale = LogNorm()
        elif scale == 'linear':
            scale = None

        pc = a.pcolormesh(Y, X, Z, vmin=vmin, vmax=vmax, norm=scale, cmap=cmap)
        a.set_yscale('linear')
        a.set_xscale('log')
        a.set_xlim((self.bins[0], self.bins[-1]))
        a.set_ylabel('Altitude (m)')
        a.set_ylim((self.layercenters[0], self.layercenters[-1]))

        a.set_xlabel('Diameter (nm)')

        a.get_yaxis().set_tick_params(direction='out', which='both')
        a.get_xaxis().set_tick_params(direction='out', which='both')

        if colorbar:
            cb = f.colorbar(pc)
            label = get_label(self.distributionType)
            cb.set_label(label)
        else:
            cb = None

        # if self.distributionType != 'calibration':
        #     a.xaxis.set_minor_formatter(plt.FormatStrFormatter("%i"))
        #     a.xaxis.set_major_formatter(plt.FormatStrFormatter("%i"))
        #
        #     f.canvas.draw()  # this is important, otherwise the ticks (at least in case of minor ticks) are not created yet
        #     ticks = a.xaxis.get_minor_ticks()
        #     for i in ticks:
        #         if i.label.get_text() in removeTickLabels:
        #             i.label.set_visible(False)

        if fit_pos:
            f, a, _ = self.normal_distribution_fits.plot_fitres(show_amplitude=False, show_width=False, ax=a)
            g = a.get_lines()[-1]
            g.set_color('m')
            g.set_label('normal dist. center')
            # a.plot(self.data.index, self.data_fit_normal.Pos, color='m', linewidth=2, label='normal dist. center')
            leg = a.legend(loc=1, fancybox=True, framealpha=0.5)
            leg.draw_frame(True)
            # if 'data_fit_normal' in dir(self):
            #     a.plot(self.data_fit_normal.Pos, self.layercenters, color='m', linewidth=2, label='normal dist. center')
            #     leg = a.legend(loc = 1, fancybox=True, framealpha=0.5)
            #     leg.draw_frame(True)

        return f, a, pc, cb

    def plot_overview(self, layers=None, show_center_of_layers = True, fit_pos = True):
        """Plot 3 plots: Size distribution Vertical profile, average size distribution, particle concentration.
         Optional layers can be defined that show up in the average plot instead of the overall average.

         Parameters
         ----------
         layers: dict, e.g. {'bottom': [0, 300]}
            define layers do average over here
         show_center_of_layers: bool
            if to show the centers of thelayers in the other plots

        """
        num_g = 10
        smallx = 3
        smally = 3

        gs = gridspec.GridSpec(num_g, num_g)
        gs.update(hspace=0.0)
        gs.update(wspace=0.0)
        f = _plt.figure()
        f.set_figwidth(f.get_figwidth() * 1.3)
        f.set_figheight(f.get_figheight() * 1.3)

        a_top = f.add_subplot(gs[:smally, :num_g - smallx])
        a_center = f.add_subplot(gs[smally:, :num_g - smallx], sharex=a_top)
        a_right = f.add_subplot(gs[smally:, num_g - smallx:], sharey=a_center)

        ####
        ## TOP plot
        if layers:
            for ln, lb in sorted(layers.items(), key=lambda x: x[1], reverse=True):
                dist_avg = self.zoom_altitude(*lb).average_overAllAltitudes().convert2dNdlogDp()
                dist_avg.plot(ax=a_top, label='{}-{}'.format(*lb))
                g = a_top.get_lines()[-1]
                col = g.get_color()
                if show_center_of_layers:
                    a_center.axhline((lb[0] + lb[1]) / 2, color = col, ls = '--')
                    a_right.axhline((lb[0] + lb[1]) / 2, color=col, ls = '--')

            leg = a_top.legend(loc = 'upper left', bbox_to_anchor=(1, 1), fontsize='x-small')
            leg.set_title('Altitude (m)')


        else:
            dist_avg = self.average_overAllAltitudes()
            dist_avg.plot(ax=a_top, label='average')

        a_top.set_yscale('log')
        # a_top.set_xlim((140, 3000))
        a_top.grid(False)

        a_top.set_yscale('log')

        ###############
        out = self.plot(ax=a_center, colorbar=False, fit_pos=fit_pos)
        pc = out[2]

        #################
        self.particle_number_concentration.plot(ax=a_right)
        a_right.set_ylabel('')
        _plt.setp(a_right.get_yticklabels(), visible=False)
        # ma_loc = MultipleLocator(20)
        # mi_loc = MultipleLocator(10)
        a_right.xaxis.set_major_locator(_MaxNLocator(4, prune='both'))
        # a_right.xaxis.set_major_locator(ma_loc)
        # a_right.xaxis.set_minor_locator(mi_loc)
        a_right.set_xlabel('Concentration (cm$^{-3}$)')
        a_right.xaxis.set_tick_params(which='major', pad=8)
        a_right.grid(False)
        a_right.set_xlim(left=5)

        cb = f.colorbar(pc)
        cb.set_label('$\\mathrm{d}N\\,/\\,\\mathrm{d}log(D_{P})$ (cm$^{-3}$)')
        return f, (a_top, a_center, a_right)


    # def plot_fitres(self, amp=True, rotate=True):
    #     """ Plots the results from fit_normal
    #
    #     Arguments
    #     ---------
    #     amp: bool.
    #         if the amplitude is to be plotted
    #     """
    #
    #     f, a = plt.subplots()
    #     a.fill_between(self.layercenters, self.data_fit_normal.Sigma_high, self.data_fit_normal.Sigma_low,
    #                    color=plt_tools.color_cycle[0],
    #                    alpha=0.5,
    #                    )
    #
    #     self.data_fit_normal.Pos.plot(ax=a, color=plt_tools.color_cycle[0], linewidth=2)
    #     g = a.get_lines()[-1]
    #     g.set_label('Center of norm. dist.')
    #     a.legend(loc=2)
    #
    #     a.set_ylabel('Particle diameter (nm)')
    #     a.set_xlabel('Altitude (m)')
    #
    #     if amp:
    #         a2 = a.twinx()
    #         self.data_fit_normal.Amp.plot(ax=a2, color=plt_tools.color_cycle[1], linewidth=2)
    #         g = a2.get_lines()[-1]
    #         g.set_label('Amplitude of norm. dist.')
    #         a2.legend()
    #         a2.set_ylabel('Amplitude - %s' % (get_label(self.distributionType)))
    #     else:
    #         a2 = False
    #     return f, a, a2

    def plot_angstromex_fit(self):
        if 'angstromexp_fit' not in dir(self):
            raise ValueError('Execute function calculate_angstromex first!')

        f, a = _plt.subplots()
        a.plot(self.angstromexp_fit.index, self.angstromexp_fit.data, 'o', color=plt_tools.color_cycle[0],
               label='exp. data')
        a.plot(self.angstromexp_fit.index, self.angstromexp_fit.fit, color=plt_tools.color_cycle[1], label='fit',
               linewidth=2)
        a.set_xlim((self.angstromexp_fit.index.min() * 0.95, self.angstromexp_fit.index.max() * 1.05))
        a.set_ylim((self.angstromexp_fit.data.min() * 0.95, self.angstromexp_fit.data.max() * 1.05))
        a.set_xlabel('Wavelength (nm)')
        a.set_ylabel('AOD')
        a.loglog()
        a.xaxis.set_minor_formatter(_plt.FormatStrFormatter("%i"))
        a.yaxis.set_minor_formatter(_plt.FormatStrFormatter("%.2f"))
        return a

    def plot_angstromex_LS(self, corr_coeff=False, std=False):
        if 'angstromexp_fit' not in dir(self):
            raise ValueError('Execute function calculate_angstromex first!')

        f, a = _plt.subplots()
        a.plot(self.angstromexp_LS.index, self.angstromexp_LS.ang_exp, color=plt_tools.color_cycle[0], linewidth=2,
               label='Angstrom exponent')
        a.set_xlabel('Altitude (m)')
        a.set_ylabel('Angstrom exponent')

        if corr_coeff:
            a.legend(loc=2)
            a2 = a.twinx()
            a2.plot(self.angstromexp_LS.index, self.angstromexp_LS.correlation_coef, color=plt_tools.color_cycle[1],
                    linewidth=2, label='corr_coeff')
            a2.set_ylabel('Correlation coefficiant')
            a2.legend(loc=1)

        if std:
            a.legend(loc=2)
            a2 = a.twinx()
            a2.plot(self.angstromexp_LS.index, self.angstromexp_LS.standard_dif, color=plt_tools.color_cycle[1],
                    linewidth=2, label='corr_coeff')
            a2.set_ylabel('Standard deviation')
            a2.legend(loc=1)

        tmp = (self.angstromexp_LS.index.max() - self.angstromexp_LS.index.min()) * 0.05
        a.set_xlim((self.angstromexp_LS.index.min() - tmp, self.angstromexp_LS.index.max() + tmp))
        return a

    def zoom_altitude(self, bottom, top):
        """'2014-11-24 16:02:30'"""
        dist = self.copy()
        dist.data = dist.data.truncate(before=bottom, after=top)
        where = _np.where(_np.logical_and(dist.layercenters <= top, dist.layercenters >= bottom))
        # dist.layercenters = dist.layercenters[where]
        dist.layerbounderies = dist.layerbounderies[where]
        if 'data_fit_normal' in dir(dist):
            dist.data_fit_normal = dist.data_fit_normal.iloc[where]
        return dist

    # dist = self.copy()
    #        dist.data = dist.data.truncate(before=start, after = end)
    #        return dist
    #


    def average_overAltitude(self, window='1S'):
        print('need fixn. Work around: generate dist_LS using a different resolution')
        self._update()
        return False

    #        window = window
    # self.data = self.data.resample(window, closed='right',label='right')
    #        if self.distributionType == 'calibration':
    #            self.data.values[_np.where(_np.isnan(self.data.values))] = 0
    #        return




    def average_overAllAltitudes(self):
        dataII = self.data.mean(axis=0)
        out = pd.DataFrame(dataII).T
        self._update()
        return SizeDist(out, self.bins, self.distributionType)


    def fit_normal(self):
        """ Fits a single normal distribution to each line in the data frame.

        Returns
        -------
        pandas DataFrame instance (also added to namespace as data_fit_normal)

        """

        super(SizeDist_LS, self).fit_normal()
        self.data_fit_normal.index = self.layercenters
        return self.data_fit_normal

#        singleHist = _np.zeros(self.data.shape[1])
#        for i in xrange(self.data.shape[1]):
#            line = self.data.values[:,i]
#            singleHist[i] = _np.average(line[~_np.isnan(line)])
#        return singleHist

def simulate_sizedistribution(diameter=[10, 2500], numberOfDiameters=100, centerOfAerosolMode=200,
                              widthOfAerosolMode=0.2, numberOfParticsInMode=1000):
    """generates a numberconcentration of an aerosol layer which has a gaussian shape when plottet in dN/log(Dp). 
    However, returned is a numberconcentrations (simply the number of particles in each bin, no normalization)
    Returns
        Number concentration (#)
        bin edges (nm)"""

    start = diameter[0]
    end = diameter[1]
    noOfD = numberOfDiameters
    centerDiameter = centerOfAerosolMode
    width = widthOfAerosolMode
    bins = _np.linspace(_np.log10(start), _np.log10(end), noOfD)
    binwidth = bins[1:] - bins[:-1]
    bincenters = (bins[1:] + bins[:-1]) / 2.
    # dNDlogDp = _plt.mlab.normpdf(bincenters, _np.log10(centerDiameter), width)# normpdf is deprecated in plt
    dNDlogDp = _sp.stats.norm.pdf(bincenters, _np.log10(centerDiameter), width)
    extraScale = 1
    scale = 1
    while 1:
        NumberConcent = dNDlogDp * binwidth * scale * extraScale
        if scale != 1:
            break
        else:
            scale = float(numberOfParticsInMode) / NumberConcent.sum()

    binEdges = 10 ** bins
    diameterBinwidth = binEdges[1:] - binEdges[:-1]

    cols = []
    for e, i in enumerate(binEdges[:-1]):
        cols.append(str(i) + '-' + str(binEdges[e + 1]))

    data = pd.DataFrame(_np.array([NumberConcent / diameterBinwidth]), columns=cols)

    return SizeDist(data, binEdges, 'dNdDp')


def simulate_sizedistribution_timeseries(diameter=[10, 2500], numberOfDiameters=100, centerOfAerosolMode=200,
                                         widthOfAerosolMode=0.2, numberOfParticsInMode=1000,
                                         startDate='2014-11-24 17:00:00',
                                         endDate='2014-11-24 18:00:00',
                                         frequency=10):
    delta = datetime.datetime.strptime(endDate, '%Y-%m-%d %H:%M:%S') - datetime.datetime.strptime(startDate,
                                                                                                  '%Y-%m-%d %H:%M:%S')
    periods = int(delta.total_seconds() / float(frequency))
    rng = pd.date_range(startDate, periods=periods, freq='%ss' % frequency)

    noOfOsz = 5
    ampOfOsz = 100

    oszi = _np.linspace(0, noOfOsz * 2 * _np.pi, periods)

    sdArray = _np.zeros((periods, numberOfDiameters - 1))
    for e, i in enumerate(rng):
        sdtmp = simulate_sizedistribution(diameter=diameter,
                                          numberOfDiameters=numberOfDiameters,
                                          widthOfAerosolMode= widthOfAerosolMode,
                                          numberOfParticsInMode = numberOfParticsInMode,
                                          centerOfAerosolMode=centerOfAerosolMode + (ampOfOsz * _np.sin(oszi[e])))
        sdArray[e] = sdtmp.data
    sdts = pd.DataFrame(sdArray, index=rng, columns=sdtmp.data.columns)
    ts = SizeDist_TS(sdts, sdtmp.bins, sdtmp.distributionType)
    ts._data_period = frequency
    return ts


def simulate_sizedistribution_layerseries(diameter=[10, 2500], numberOfDiameters=100, heightlimits=[0, 6000],
                                          noOflayers=100, layerHeight=[500., 4000.], layerThickness=[100., 300.],
                                          layerDensity=[1000., 5000.], layerModecenter=[200., 800.], widthOfAerosolMode = 0.2 ):
    gaussian = lambda x, mu, sig: _np.exp(-(x - mu) ** 2 / (2 * sig ** 2))

    lbt = _np.linspace(heightlimits[0], heightlimits[1], noOflayers + 1)
    layerbounderies = _np.array([lbt[:-1], lbt[1:]]).transpose()
    layercenter = (lbt[1:] + lbt[:-1]) / 2.

    # strata = _np.linspace(heightlimits[0],heightlimits[1],noOflayers+1)

    layerArray = _np.zeros((noOflayers, numberOfDiameters - 1))

    for e, stra in enumerate(layercenter):
        for i, lay in enumerate(layerHeight):
            sdtmp = simulate_sizedistribution(diameter=diameter, numberOfDiameters=numberOfDiameters,
                                              widthOfAerosolMode=widthOfAerosolMode, centerOfAerosolMode=layerModecenter[i],
                                              numberOfParticsInMode=layerDensity[i])
            layerArray[e] += sdtmp.data.values[0] * gaussian(stra, layerHeight[i], layerThickness[i])

    sdls = pd.DataFrame(layerArray, index=layercenter, columns=sdtmp.data.columns)
    return SizeDist_LS(sdls, sdtmp.bins, sdtmp.distributionType, layerbounderies)


def generate_aerosolLayer(diameter=[.01, 2.5], numberOfDiameters=30, centerOfAerosolMode=0.6,
                          widthOfAerosolMode=0.2, numberOfParticsInMode=10000, layerBoundery=[0., 10000], ):
    """Probably deprecated!?! generates a numberconcentration of an aerosol layer which has a gaussian shape when plottet in dN/log(Dp). 
    However, returned is a numberconcentrations (simply the number of particles in each bin, no normalization)
    Returns
        Number concentration (#)
        bin edges (nm)"""

    layerBoundery = _np.array(layerBoundery)
    start = diameter[0]
    end = diameter[1]
    noOfD = numberOfDiameters
    centerDiameter = centerOfAerosolMode
    width = widthOfAerosolMode
    bins = _np.linspace(_np.log10(start), _np.log10(end), noOfD)
    binwidth = bins[1:] - bins[:-1]
    bincenters = (bins[1:] + bins[:-1]) / 2.
    dNDlogDp = _plt.mlab.normpdf(bincenters, _np.log10(centerDiameter), width)
    extraScale = 1
    scale = 1
    while 1:
        NumberConcent = dNDlogDp * binwidth * scale * extraScale
        if scale != 1:
            break
        else:
            scale = float(numberOfParticsInMode) / NumberConcent.sum()

    binEdges = 10 ** bins
    # diameterBinCenters = (binEdges[1:] + binEdges[:-1])/2.
    diameterBinwidth = binEdges[1:] - binEdges[:-1]

    cols = []
    for e, i in enumerate(binEdges[:-1]):
        cols.append(str(i) + '-' + str(binEdges[e + 1]))

    layerBoundery = _np.array([0., 10000.])
    # layerThickness = layerBoundery[1:] - layerBoundery[:-1]
    layerCenter = [5000.]
    data = pd.DataFrame(_np.array([NumberConcent / diameterBinwidth]), index=layerCenter, columns=cols)
    # return data

    #     atmosAerosolNumberConcentration = pd.DataFrame()
    # atmosAerosolNumberConcentration['bin_center'] = pd.Series(diameterBinCenters)
    #     atmosAerosolNumberConcentration['bin_start'] = pd.Series(binEdges[:-1])
    #     atmosAerosolNumberConcentration['bin_end'] = pd.Series(binEdges[1:])
    #     atmosAerosolNumberConcentration['numberConcentration'] = pd.Series(NumberConcent)
    #     return atmosAerosolNumberConcentration

    return SizeDist_LS(data, binEdges, 'dNdDp', layerBoundery)


def test_generate_numberConcentration():
    """result should look identical to Atmospheric Chemistry and Physis page 422"""
    nc = generate_aerosolLayer(diameter=[0.01, 10], centerOfAerosolMode=0.8, widthOfAerosolMode=0.3,
                               numberOfDiameters=100, numberOfParticsInMode=1000, layerBoundery=[0.0, 10000])

    _plt.plot(nc.bincenters, nc.data.values[0].transpose() * nc.binwidth, label='numberConc')
    _plt.plot(nc.bincenters, nc.data.values[0].transpose(), label='numberDist')
    ncLN = nc.convert2dNdlogDp()
    _plt.plot(ncLN.bincenters, ncLN.data.values[0].transpose(), label='LogNormal')
    _plt.legend()
    _plt.semilogx()








def test_ext_coeff_vertical_profile():
    #todo: make this a real test
    dist = simulate_sizedistribution_layerseries(layerHeight=[3000.0, 3000.0],
                                               layerDensity=[1000.0, 100.0],
                                               layerModecenter=[100.0, 100.0],
                                               layerThickness=[6000, 6000],
                                               widthOfAerosolMode = 0.01,
                                               noOflayers=3,
                                               numberOfDiameters=1000)
    dist.plot()


    dist = dist.zoom_diameter(99,101)
    avg = dist.average_overAllAltitudes()
    f,a = avg.plot()
    a.set_xscale('linear')

    opt = dist.calculate_optical_properties(550, n = 1.455)
    opt_II = dist.calculate_optical_properties(550, n = 1.1)
    opt_III = dist.calculate_optical_properties(550, n = 4.)

    ext = opt.get_extinction_coeff_verticle_profile()
    ext_II = opt_II.get_extinction_coeff_verticle_profile()
    ext_III = opt_III.get_extinction_coeff_verticle_profile()

    tvI_is = (ext_III.data/ext.data).values[0][0]
    tvI_want = 14.3980239083
    tvII_is = (ext_II.data/ext.data).values[0][0]
    tvII_want = 0.05272993413

    print('small deviations could come from averaging over multiple bins with slightly different diameter')
    print('test values 1 is/should_be: %s/%s'%(tvI_is,tvI_want))
    print('test values 2 is/should_be: %s/%s'%(tvII_is,tvII_want))

    return False