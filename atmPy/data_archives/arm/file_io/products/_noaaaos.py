import warnings as _warnings
from atmPy.general import timeseries as _timeseries
# from atmPy.data_archives.arm import _tools
import pandas as _pd
from atmPy.data_archives.arm._netCDF import ArmDataset as _ArmDataset
import numpy as _np
from atmPy.tools import decorators
from atmPy.aerosols.physics import hygroscopicity as hygrow
from atmPy.aerosols.physics import hygroscopicity
import atmPy.aerosols.instruments.nephelometer.nephelometer as nephelometer
import atmPy.aerosols.instruments.photometer.photometer as photometer
try:
    import xarray as _xr
except ModuleNotFoundError:
    _warnings.warn('netCDF4 not installed. You might encounter some functionality limitations.')
import os as _os
from atmPy.data_archives.arm import _tools
# data_archives.arm._tools as _tools

# import pdb as _pdb


def read_noaaaos(fname, in_time_window = None, keep_xr_dataset = False, verbose = False):
    noaaaos = Noaaaos()

    noaaaos.nephelometer_1um = nephelometer.Nephelometer()
    noaaaos.nephelometer_10um = nephelometer.Nephelometer()
    noaaaos.psap_1um = photometer.Photometer()
    noaaaos.psap_10um = photometer.Photometer()


    def get_scatt_coeff(ds):
        df = _pd.DataFrame()
        df['450 nm'] = ds.Bs_B_Dry_1um_Neph3W_1.to_pandas()
        df['550 nm'] = ds.Bs_G_Dry_1um_Neph3W_1.to_pandas()
        df['700 nm'] = ds.Bs_R_Dry_1um_Neph3W_1.to_pandas()
        df[df == -9999.0] = _np.nan
        scatt_coeff_1um_dry = df

        df = _pd.DataFrame()
        df['450 nm'] = ds.Bs_B_Dry_10um_Neph3W_1.to_pandas()
        df['550 nm'] = ds.Bs_G_Dry_10um_Neph3W_1.to_pandas()
        df['700 nm'] = ds.Bs_R_Dry_10um_Neph3W_1.to_pandas()
        df[df == -9999.0] = _np.nan
        scatt_coeff_10um_dry = df

        return scatt_coeff_1um_dry, scatt_coeff_10um_dry

    def get_scatt_coeff_variability(ds):
        df = _pd.DataFrame()
        df['450 nm'] = ds.Bs_B_1um_PSAP1W_1_std.to_pandas()
        df['550 nm'] = ds.Bs_G_1um_PSAP1W_1_std.to_pandas()
        df['700 nm'] = ds.Bs_R_1um_PSAP1W_1_std.to_pandas()
        df[df == -9999.0] = _np.nan
        scatt_coeff_1um_dry = df

        df = _pd.DataFrame()
        df['450 nm'] = ds.Bs_B_Dry_10um_PSAP1W_1_std.to_pandas()
        df['550 nm'] = ds.Bs_G_Dry_10um_PSAP1W_1_std.to_pandas()
        df['700 nm'] = ds.Bs_R_Dry_10um_PSAP1W_1_std.to_pandas()
        df[df == -9999.0] = _np.nan
        scatt_coeff_10um_dry = df
        return scatt_coeff_1um_dry, scatt_coeff_10um_dry

    def get_hembackscatt(ds):
        df = _pd.DataFrame()
        df['450 nm'] = ds.Bbs_B_Dry_1um_Neph3W_1.to_pandas()
        df['550 nm'] = ds.Bbs_G_Dry_1um_Neph3W_1.to_pandas()
        df['700 nm'] = ds.Bbs_R_Dry_1um_Neph3W_1.to_pandas()
        df[df == -9999.0] = _np.nan
        scatt_coeff_1um_dry = df

        df = _pd.DataFrame()
        df['450 nm'] = ds.Bbs_B_Dry_10um_Neph3W_1.to_pandas()
        df['550 nm'] = ds.Bbs_G_Dry_10um_Neph3W_1.to_pandas()
        df['700 nm'] = ds.Bbs_R_Dry_10um_Neph3W_1.to_pandas()
        df[df == -9999.0] = _np.nan
        scatt_coeff_10um_dry = df
        return scatt_coeff_1um_dry, scatt_coeff_10um_dry

    if type(fname) == str:
        if _os.path.isdir(fname):
            # fname = fname + _os.listdir(fname)
            fname = [fname + i for i in _os.listdir(fname)]
        elif _os.path.isfile(fname):
            fname = [fname]

    neph_scatt_coeff_1um = []
    neph_scatt_coeff_10um = []
    neph_scatt_coeff_1um_std = []
    neph_scatt_coeff_10um_std = []
    neph_back_scatt_coeff_1um = []
    neph_back_scatt_coeff_10um = []
    neph_RH = []
    psap_abs_1um = []
    psap_abs_10um = []

    found_one_in_window = False
    for fn in fname:
        if verbose:
            print('reading: {}'.format(fn))
        if in_time_window:
            if not _tools.is_in_time_window(fn, in_time_window, verbose = verbose):
                continue
            else:
                found_one_in_window = True
        ds = _xr.open_dataset(fn)
        version = ds.attrs['zeb_platform']

        sc1, sc10 = get_scatt_coeff(ds)
        neph_scatt_coeff_1um.append(sc1)
        neph_scatt_coeff_10um.append(sc10)

        if version in []:
            sc1, sc10 = get_scatt_coeff_variability(ds)
            neph_scatt_coeff_1um_std.append(sc1)
            neph_scatt_coeff_10um_std.append(sc10)

        sc1, sc10 = get_hembackscatt(ds)
        neph_back_scatt_coeff_1um.append(sc1)
        neph_back_scatt_coeff_10um.append(sc10)

        df = _pd.DataFrame()
        df['565 nm'] = ds.Ba_G_Dry_1um_PSAP1W_1.to_pandas()
        psap_abs_1um.append(df)

        df = _pd.DataFrame()
        df['565 nm'] = ds.Ba_G_Dry_10um_PSAP1W_1.to_pandas()
        psap_abs_10um.append(df)

    if in_time_window and not found_one_in_window:
        raise KeyError('No file found in defind time window (in_time_window = {})'.format(in_time_window))

    noaaaos.nephelometer_1um.scattering_coeff = _timeseries.TimeSeries(_pd.concat(neph_scatt_coeff_1um).sort_index())
    noaaaos.nephelometer_10um.scattering_coeff = _timeseries.TimeSeries(_pd.concat(neph_scatt_coeff_10um).sort_index())

    if version in []:
        noaaaos.nephelometer_1um.scattering_coeff.standard_diviation = _timeseries.TimeSeries(_pd.concat(neph_scatt_coeff_1um_std).sort_index())
        noaaaos.nephelometer_10um.scattering_coeff.standard_diviation = _timeseries.TimeSeries(_pd.concat(neph_scatt_coeff_10um_std).sort_index())

    noaaaos.nephelometer_1um.hemisphericbackscatt_coeff = _timeseries.TimeSeries(_pd.concat(neph_back_scatt_coeff_1um).sort_index())
    noaaaos.nephelometer_10um.hemisphericbackscatt_coeff = _timeseries.TimeSeries(_pd.concat(neph_back_scatt_coeff_10um).sort_index())

    noaaaos.psap_1um.absorbtion_coeff = _timeseries.TimeSeries(_pd.concat(psap_abs_1um).sort_index())
    noaaaos.psap_10um.absorbtion_coeff = _timeseries.TimeSeries(_pd.concat(psap_abs_10um).sort_index())
    if keep_xr_dataset:
        noaaaos.xr_dataset = ds
    return noaaaos


class Noaaaos(object):
    def __init__(self):
        self._nephelometer_1um = None
        self._nephelometer_10um = None
        self._psap_1um = None
        self._psap_10um = None
        self._cpc = None

    @property
    def psap_1um(self):
        return self._psap_1um

    @psap_1um.setter
    def psap_1um(self, value):
        self._psap_1um = value

    @property
    def psap_10um(self):
        return self._psap_10um

    @psap_10um.setter
    def psap_10um(self, value):
        self._psap_10um = value

    @property
    def cpc(self):
        return self._cpc

    @cpc.setter
    def cpc(self, value):
        self._cpc = value

#######################################################################################
######## Older stuff
#######################################################################################
def calculate_f_RH(noaaaos, RH_center, RH_tolerance, which):
    """

    Parameters
    ----------
    noaaaos: noaaaos.ArmDataset instance
    RH_center: int
        The wet nephelometer RH value varies. This is the RH value at which the
        neph data is used. (typical value is 85)
    RH_tolerance: float
        Defines the range of data around RH_center that is included,
        range = RH_center +- RH_tolerance
        Typical value for RH_tolerance is 1
    which: str
        The nephelometer has 3 wavelength channels. Choose between:
        "all", "gree", "red", or "blue".


    Returns
    -------
    TimeSeries instance

    """
    lim = (RH_center - RH_tolerance, RH_center + RH_tolerance)
    rh_wet = noaaaos.RH_nephelometer.data.RH_NephVol_Wet.copy()
    rh_wet.values[rh_wet.values > lim[1]] = _np.nan
    rh_wet.values[rh_wet.values < lim[0]] = _np.nan
# rh_wet.plot()
    if which == 'all':
        which = ['green', 'red', 'blue']
    else:
        which = [which]

    df = _pd.DataFrame(index = noaaaos.scatt_coeff.data.index)
    for col in which:
        if col == 'green':
            wet_column = 'Bs_G_Wet_1um_Neph3W_2'
            dry_column = 'Bs_G_Dry_1um_Neph3W_1'
        elif col == 'red':
            wet_column = 'Bs_R_Wet_1um_Neph3W_2'
            dry_column = 'Bs_R_Dry_1um_Neph3W_1'
        elif col == 'blue':
            wet_column = 'Bs_B_Wet_1um_Neph3W_2'
            dry_column = 'Bs_B_Dry_1um_Neph3W_1'
        else:
            txt = '%s is not an option. Choose between ["all", "green", "red", "blue"]'
            raise ValueError(txt)


        scatt_coeff_wet = noaaaos.scatt_coeff.data[wet_column].copy()
        scatt_coeff_dry = noaaaos.scatt_coeff.data[dry_column].copy()

        scatt_coeff_wet.values[_np.isnan(rh_wet.values)] = _np.nan
        scatt_coeff_dry.values[_np.isnan(scatt_coeff_wet.values)] = _np.nan
        f_rh = scatt_coeff_wet/scatt_coeff_dry

        f_rh_int = f_rh.interpolate()
        f_rh_mean = _pd.rolling_mean(f_rh_int, 40, center = True)
        df[col] = f_rh_mean
    ts = _timeseries.TimeSeries(df)
    ts._y_label = '$f(RH = %i \pm %i \%%)$'%(RH_center, RH_tolerance)
    return ts




class ArmDatasetSub(_ArmDataset):
    def __init__(self,*args, **kwargs):
        self._data_period = 60
        self._time_offset = (- self._data_period, 's')
        super(ArmDatasetSub,self).__init__(*args, **kwargs)

        self._concatable = ['abs_coeff', 'back_scatt', 'scatt_coeff', 'RH_nephelometer']

        self.__f_of_RH = None
        self.__kappa = None
        self.__growthfactor = None
        self.__hemispheric_backscattering_ratio = None

        self.__sup_fofRH_RH_center = None
        self.__sup_fofRH_RH_tolerance = None
        self.__sup_fofRH_which = None
        self.__sup_kappa_sizedist = None
        self.__sup_kappa_wavelength = None
        self.__hygroscopicity = None
        self.__hygroscopicity_red = None
        self.__hygroscopicity_blue = None



    def _parse_netCDF(self):

        abs_coeff = ['Ba_G_Dry_10um_PSAP1W_1',
                    'Ba_G_Dry_1um_PSAP1W_1',
                    'Ba_B_Dry_10um_PSAP3W_1',
                    'Ba_G_Dry_10um_PSAP3W_1',
                    'Ba_R_Dry_10um_PSAP3W_1',
                    'Ba_B_Dry_1um_PSAP3W_1',
                    'Ba_G_Dry_1um_PSAP3W_1',
                    'Ba_R_Dry_1um_PSAP3W_1',
                    ]

        scat_coeff =   ['Bs_B_Dry_10um_Neph3W_1',
                            'Bs_G_Dry_10um_Neph3W_1',
                            'Bs_R_Dry_10um_Neph3W_1',
                            'Bs_B_Wet_10um_Neph3W_2',
                            'Bs_G_Wet_10um_Neph3W_2',
                            'Bs_R_Wet_10um_Neph3W_2',
                            'Bs_B_Dry_1um_Neph3W_1',
                            'Bs_G_Dry_1um_Neph3W_1',
                            'Bs_R_Dry_1um_Neph3W_1',
                            'Bs_B_Wet_1um_Neph3W_2',
                            'Bs_G_Wet_1um_Neph3W_2',
                            'Bs_R_Wet_1um_Neph3W_2',
                            ]


        bscat_coeff_vars = ['Bbs_B_Dry_10um_Neph3W_1',
                            'Bbs_G_Dry_10um_Neph3W_1',
                            'Bbs_R_Dry_10um_Neph3W_1',
                            'Bbs_B_Wet_10um_Neph3W_2',
                            'Bbs_G_Wet_10um_Neph3W_2',
                            'Bbs_R_Wet_10um_Neph3W_2',
                            'Bbs_B_Dry_1um_Neph3W_1',
                            'Bbs_G_Dry_1um_Neph3W_1',
                            'Bbs_R_Dry_1um_Neph3W_1',
                            'Bbs_B_Wet_1um_Neph3W_2',
                            'Bbs_G_Wet_1um_Neph3W_2',
                            'Bbs_R_Wet_1um_Neph3W_2',
                            ]

        RH_neph = ['RH_NephVol_Dry',
              'RH_NephVol_Wet']

        def var2ts(self, var_list, column_name):
            """extracts the list of variables from the file_obj and puts them all in one data frame"""
            df = _pd.DataFrame(index = self.time_stamps)
            for var in var_list:
                data = self._read_variable(var)
                df[var] = _pd.Series(data['data'], index = self.time_stamps)
            df.columns.name = column_name
            out = _timeseries.TimeSeries(df)
            out._data_period = self._data_period
            return out
        self.abs_coeff = var2ts(self, abs_coeff, 'abs_coeff_1/Mm')
        self.scatt_coeff = var2ts(self, scat_coeff, 'scatt_coeff_1/Mm')
        self.back_scatt = var2ts(self, bscat_coeff_vars, 'back_scatt_1/Mm')
        self.RH_nephelometer = var2ts(self, RH_neph, 'RH')


    def plot_all(self):
        self.abs_coeff.plot()
        self.back_scatt.plot()
        self.scatt_coeff.plot()
        self.RH_nephelometer.plot()

    #########
    ### OLD - use hygroscopicity
    # @property
    # @decorators.change_doc(calculate_f_RH)
    # def f_of_RH(self):
    #     """
    #     Parameter names are changed to self.sup_fofRH_RH_center, self.sup_fofRH_RH_tolerance, and self.sup_fofRH_which
    #     """
    #     if not self.__f_of_RH:
    #         if not self.sup_fofRH_RH_center or not self.sup_fofRH_RH_tolerance or not self.sup_fofRH_which:
    #             txt = "Make sure you define the following attributes first: \nself.sup_fofRH_RH_center, self.sup_fofRH_RH_tolerance, self.sup_fofRH_which"
    #             raise ValueError(txt)
    #         self.__f_of_RH = calculate_f_RH(self,self.sup_fofRH_RH_center, self.sup_fofRH_RH_tolerance, self.sup_fofRH_which)
    #         self.__f_of_RH._data_period = self._data_period
    #     return self.__f_of_RH
    #
    # @f_of_RH.setter
    # def f_of_RH(self, value):
    #     self.__f_of_RH = value


    @property
    def hemispheric_backscattering_ratio(self):
        if not self.__hemispheric_backscattering_ratio:
            if _np.any(self.back_scatt.data.index != self.scatt_coeff.data.index):
                raise IndexError(
                    "The indeces doe not seam to match, that should not be possible!")

            bdf = self.back_scatt.data.copy()
            sdf = self.scatt_coeff.data.copy()

            bk = [i.replace('Bbs_', '') for i in bdf.keys()]
            sk = [i.replace('Bs_', '') for i in sdf.keys()]
            if bk != sk:
                raise KeyError(
                    'These two data frames seam to be not the right ones ... headers do not match (%s,%s)' % (
                    bk, sk))

            new_col_names = bk
            bdf.columns = new_col_names
            sdf.columns = new_col_names

            out = _timeseries.TimeSeries(bdf.div(sdf))
            out._data_period = self.back_scatt._data_period
            self.__hemispheric_backscattering_ratio = out

        return self.__hemispheric_backscattering_ratio

    @hemispheric_backscattering_ratio.setter
    def hemispheric_backscattering_ratio(self,value):
        self.__hemispheric_backscattering_ratio = value

    @property
    @decorators.change_doc(hygrow.kappa_from_fofrh_and_sizedist)
    def kappa(self):
        if not self.__kappa:
            if not self.sup_kappa_sizedist or not self.sup_kappa_wavelength:
                txt = "Make sure you define the following attributes first: \nself.sup_kappa_sizedist and self.sup_kappa_wavelength"
                raise ValueError(txt)

            self.__kappa, self.__growthfactor = hygrow.kappa_from_fofrh_and_sizedist(self.f_of_RH, self.sup_kappa_sizedist,
                                                                                     self.sup_kappa_wavelength, self.sup_fofRH_RH_center)
        return self.__kappa

    @kappa.setter
    def kappa(self,value):
        self.__kappa = value

    @property
    @decorators.change_doc(hygrow.kappa_from_fofrh_and_sizedist)
    def growth_factor(self):
        if self.__growthfactor:
            self.kappa
        return self.__growthfactor

    @property
    def sup_kappa_wavelength(self):
        return self.__sup_kappa_wavelength

    @sup_kappa_wavelength.setter
    def sup_kappa_wavelength(self, value):
        self.__kappa = None
        self.__growthfactor = None
        self.__sup_kappa_wavelength = value

    @property
    def sup_kappa_sizedist(self):
        return self.__sup_kappa_sizedist

    @sup_kappa_sizedist.setter
    def sup_kappa_sizedist(self, value):
        self.__kappa = None
        self.__growthfactor = None
        self.__sup_kappa_sizedist = value


    @property
    def hygroscopicity_10um(self):
        if not self.__hygroscopicity:
            self.__hygroscopicity = hygroscopicity.fofRH_from_dry_wet_scattering(self.scatt_coeff._del_all_columns_but('Bs_G_Dry_10um_Neph3W_1'),
                                                                                 self.scatt_coeff._del_all_columns_but('Bs_G_Wet_10um_Neph3W_2'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Dry'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Wet'),
                                                                                 return_fits = False)
        return self.__hygroscopicity

    @property
    def hygroscopicity_10um_red(self):
        if not self.__hygroscopicity_red:
            self.__hygroscopicity_red = hygroscopicity.fofRH_from_dry_wet_scattering(self.scatt_coeff._del_all_columns_but('Bs_R_Dry_10um_Neph3W_1'),
                                                                                 self.scatt_coeff._del_all_columns_but('Bs_R_Wet_10um_Neph3W_2'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Dry'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Wet'),
                                                                                 return_fits = False)
        return self.__hygroscopicity_red

    @property
    def hygroscopicity_10um_blue(self):
        if not self.__hygroscopicity_blue:
            self.__hygroscopicity_blue = hygroscopicity.fofRH_from_dry_wet_scattering(self.scatt_coeff._del_all_columns_but('Bs_B_Dry_10um_Neph3W_1'),
                                                                                 self.scatt_coeff._del_all_columns_but('Bs_B_Wet_10um_Neph3W_2'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Dry'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Wet'),
                                                                                 return_fits = False)
        return self.__hygroscopicity_blue

    @property
    def hygroscopicity_1um(self):
        if not self.__hygroscopicity:
            self.__hygroscopicity = hygroscopicity.fofRH_from_dry_wet_scattering(self.scatt_coeff._del_all_columns_but('Bs_G_Dry_1um_Neph3W_1'),
                                                                                 self.scatt_coeff._del_all_columns_but('Bs_G_Wet_1um_Neph3W_2'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Dry'),
                                                                                 self.RH_nephelometer._del_all_columns_but('RH_NephVol_Wet'),
                                                                                 return_fits=False)
        return self.__hygroscopicity


    # @property
    # def sup_fofRH_which(self):
    #     return self.__sup_fofRH_which
    #
    # @sup_fofRH_which.setter
    # def sup_fofRH_which(self, value):
    #     self.__f_of_RH = None
    #     self.__sup_fofRH_which = value
    #
    # @property
    # def sup_fofRH_RH_center(self):
    #     return self.__sup_fofRH_RH_center
    #
    # @sup_fofRH_RH_center.setter
    # def sup_fofRH_RH_center(self, value):
    #     self.__f_of_RH = None
    #     self.__sup_fofRH_RH_center = value
    #
    # @property
    # def sup_fofRH_RH_tolerance(self):
    #     return self.__sup_fofRH_RH_tolerance
    #
    # @sup_fofRH_RH_tolerance.setter
    # def sup_fofRH_RH_tolerance(self, value):
    #     self.__f_of_RH = None
    #     self.__sup_fofRH_RH_tolerance = value


def _concat_rules(arm_data_objs):
    """nothing here"""
    # out = arm_data_obj
    out = ArmDatasetSub(False)
    out._concat(arm_data_objs)
    return out
#
# def _concat_rules(arm_data_objs):
#     out = ArmDatasetSub(False)
#     out.abs_coeff = _timeseries.TimeSeries(_pd.concat([i.abs_coeff.data for i in arm_data_objs]))
#     out.abs_coeff._data_period = out._data_period
#     out.back_scatt = _timeseries.TimeSeries(_pd.concat([i.back_scatt.data for i in arm_data_objs]))
#     out.back_scatt._data_period = out._data_period
#     out.scatt_coeff = _timeseries.TimeSeries(_pd.concat([i.scatt_coeff.data for i in arm_data_objs]))
#     out.scatt_coeff._data_period = out._data_period
#     out.RH_nephelometer = _timeseries.TimeSeries(_pd.concat([i.RH_nephelometer.data for i in arm_data_objs]))
#     out.RH_nephelometer._data_period = out._data_period
#     out.time_stamps = out.abs_coeff.data.index
#     return out
#

