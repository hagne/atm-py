from atmPy.general import timeseries as _timeseries
from atmPy.data_archives.arm import _netCDF
import pandas as _pd
import numpy as _np
from atmPy.aerosols.physics import hygroscopic_growth as _hygrow
from atmPy.tools import decorators as _decorators

class ArmDatasetSub(_netCDF.ArmDataset):

    info = ("This data product has a few gotchas:\n"
            "- The function profided is not agains 40 as it suggests, it must be some other value, which I don't what it is")

    def __init__(self,*args, **kwargs):
        self._data_period = 3600
        self._time_offset = (- self._data_period, 's')
        super(ArmDatasetSub,self).__init__(*args, **kwargs)

        self._concatable = ['f_RH_scatt_funcs_2p', 'f_RH_scatt_funcs_3p','f_RH_scatt_2p_85_40', 'f_RH_scatt_3p_85_40', 'f_RH_scatt_2p_ab_G_1um']


        ####
        # for properties
        self.__f_RH_scatt_2p = None
        self.__f_RH_scatt_3p  = None
        self.__sup_RH = None
        self.__kappa = None
        self.__growthfactor = None
        self.__kappa_85_40 = None
        self.__growthfactor_85_40 = None
        self.__sup_kappa_sizedist = None
        self.__sup_kappa_wavelength = None


    def _data_quality_control(self):
        #######
        ## Define what is good, patchy or bad data
        ## delete if not quality flags exist
        ## settings are ignored if no quality flags are profided (e.g. noaaaos)

        if self.data_quality_flag_max == None:
            if self.data_quality == 'good':
                self.data_quality_flag_max = 0
            elif self.data_quality == 'patchy':
                self.data_quality_flag_max = 7
            elif self.data_quality == 'bad':
                self.data_quality_flag_max = 100000
            else:
                txt = '%s is not an excepted values for data_quality ("good", "patchy", "bad")'%(self.data_quality)
                raise ValueError(txt)

    def _parse_netCDF(self):
        super(ArmDatasetSub,self)._parse_netCDF()
        # self.rh = self._read_variable2timeseries(['rh_60m', 'rh_60m'], column_name='Relative Humidity (%)')


        # for the 2 parameter function
        def ab_2_f_RH_func(ab):
            ab = ab.copy()
            a, b = ab
            # a = 1. # I was just told that a is supposed to be set to one from Ann (upstairs)
            f_RH = lambda RH: a * (1 - (RH / 100.)) ** (-b)  # 'bsp(RH%)/Bsp(~40%) = a*[1-(RH%/100)]^(-b)'
            return f_RH

        varies = ['fRH_Bs_R_10um_2p',
                  'fRH_Bs_G_10um_2p',
                  'fRH_Bs_B_10um_2p',
                  'fRH_Bs_R_1um_2p',
                  'fRH_Bs_G_1um_2p',
                  'fRH_Bs_B_1um_2p']

        df = _pd.DataFrame(index=self.time_stamps)
        df_ab = _pd.DataFrame(index=self.time_stamps)
        for key in varies:
            data = self._read_variable(key, reverse_qc_flag=8)
            dft = _pd.DataFrame(data, index=self.time_stamps)
            df[key] = dft.apply(ab_2_f_RH_func, axis=1)
            if key == 'fRH_Bs_G_1um_2p':
                self.f_RH_scatt_2p_ab_G_1um = _timeseries.TimeSeries(_pd.DataFrame(dft))
                self.f_RH_scatt_2p_ab_G_1um._data_period = self._data_period

        self.f_RH_scatt_funcs_2p = _timeseries.TimeSeries(df)
        self.f_RH_scatt_funcs_2p._data_period = self._data_period




        #for the 3 parameter function
        def abc_2_f_RH_func(abc):
            abc = abc.copy()
            a, b, c = abc
            # a = 1.
            f_RH = lambda RH: a * (1 + (b * (RH / 100.)**c))
            return f_RH

        varies = ['fRH_Bs_R_10um_3p',
                  'fRH_Bs_G_10um_3p',
                  'fRH_Bs_B_10um_3p',
                  'fRH_Bs_R_1um_3p',
                  'fRH_Bs_G_1um_3p',
                  'fRH_Bs_B_1um_3p']

        df = _pd.DataFrame(index=self.time_stamps)
        for key in varies:
            data = self._read_variable(key, reverse_qc_flag=8)
            dft = _pd.DataFrame(data, index=self.time_stamps)
            df[key] = dft.apply(abc_2_f_RH_func, axis=1)
        self.f_RH_scatt_funcs_3p = _timeseries.TimeSeries(df)
        self.f_RH_scatt_funcs_3p._data_period = self._data_period

        # f or RH at predifined point
        varies = ['ratio_85by40_Bs_R_10um_2p',
                  'ratio_85by40_Bs_G_10um_2p',
                  'ratio_85by40_Bs_B_10um_2p',
                  'ratio_85by40_Bs_R_1um_2p',
                  'ratio_85by40_Bs_G_1um_2p',
                  'ratio_85by40_Bs_B_1um_2p']

        self.f_RH_scatt_2p_85_40 = self._read_variable2timeseries(varies,reverse_qc_flag=8)

        varies = ['ratio_85by40_Bs_R_10um_3p',
                  'ratio_85by40_Bs_G_10um_3p',
                  'ratio_85by40_Bs_B_10um_3p',
                  'ratio_85by40_Bs_R_1um_3p',
                  'ratio_85by40_Bs_G_1um_3p',
                  'ratio_85by40_Bs_B_1um_3p']

        self.f_RH_scatt_3p_85_40 = self._read_variable2timeseries(varies, reverse_qc_flag=8)

        varies = ['ratio_85by40_Bbs_R_10um_2p',
                  'ratio_85by40_Bbs_G_10um_2p',
                  'ratio_85by40_Bbs_B_10um_2p',
                  'ratio_85by40_Bbs_R_1um_2p',
                  'ratio_85by40_Bbs_G_1um_2p',
                  'ratio_85by40_Bbs_B_1um_2p']

        self.f_RH_backscatt_2p_85_40 = self._read_variable2timeseries(varies, reverse_qc_flag=8)


    def plot_all(self):
        self.rh.plot()

    @property
    def f_RH_scatt_3p(self):
        """Note, when calculating a f(RH) with this function it has a mysterious off set in it.
        When you plan is to calculate f(RH) between 80 and 40 you actually have to apply this function for
        both values and than divide."""
        if not self.__f_RH_scatt_3p:
            if not self.sup_RH:
                raise ValueError('please set the relative humidity in sup_RH')

            def applyfunk(value):
                if type(value).__name__ == 'function':
                    return value(self.sup_RH)
                else:
                    return _np.nan

            # data = self.f_RH_scatt_funcs.data.applymap(lambda x: x(self.sup_RH))
            data = self.f_RH_scatt_funcs_3p.data.applymap(applyfunk)
            # data = _pd.DataFrame(data, columns=['f_%i'%(self.sup_RH)])
            self.__f_RH_scatt_3p = _timeseries.TimeSeries(data)
            self.__f_RH_scatt_3p._data_period = self.f_RH_scatt_funcs_3p._data_period
        return self.__f_RH_scatt_3p

    @property
    def f_RH_scatt_2p(self):
        """Note, when calculating a f(RH) with this function it has a mysterious off set in it.
        When you plan is to calculate f(RH) between 80 and 40 you actually have to apply this function for
        both values and than divide."""
        if not self.__f_RH_scatt_2p:
            if not self.sup_RH:
                raise ValueError('Please set the relative humidity to calculate f(RH) at. Therefore, set sup_RH (in %)')

            def applyfunk(value):
                if type(value).__name__ == 'function':
                    return value(self.sup_RH)
                else:
                    return _np.nan

            # data = self.f_RH_scatt_funcs.data.applymap(lambda x: x(self.sup_RH))
            data = self.f_RH_scatt_funcs_2p.data.applymap(applyfunk)
            # data = _pd.DataFrame(data, columns=['f_%i'%(self.sup_RH)])
            self.__f_RH_scatt_2p = _timeseries.TimeSeries(data)
            self.__f_RH_scatt_2p._data_period = self.f_RH_scatt_funcs_2p._data_period
        return self.__f_RH_scatt_2p

    @property
    def sup_RH(self):
        return self.__sup_RH

    @sup_RH.setter
    def sup_RH(self, value):
        self.__sup_RH = value
        self.__f_RH_scatt_3p = None
        self.__f_RH_scatt_2p = None

    @property
    @_decorators.change_doc(_hygrow.kappa_from_fofrh_and_sizedist)
    def kappa(self):
        if not self.__kappa:
            if not self.sup_kappa_sizedist or not self.sup_kappa_wavelength or not self.sup_RH:
                txt = ("Make sure you define the following attributes first: \n"
                       "\t - self.sub_RH -- needed to calculate f(RH) for a representitive value (e.g. 80)\n"
                       "\t - self.sup_kappa_sizedist\n"
                       "\t - and self.sup_kappa_wavelength")
                raise ValueError(txt)

            self.__kappa, self.__growthfactor = _hygrow.kappa_from_fofrh_and_sizedist(self.f_RH_scatt_2p, self.sup_kappa_sizedist,
                                                                                      self.sup_kappa_wavelength, self.sup_RH)
        return self.__kappa

    @kappa.setter
    def kappa(self, value):
        self.__kappa = value

    @property
    @_decorators.change_doc(_hygrow.kappa_from_fofrh_and_sizedist)
    def growth_factor(self):
        if self.__growthfactor:
            self.kappa
        return self.__growthfactor



    @property
    @_decorators.change_doc(_hygrow.kappa_from_fofrh_and_sizedist)
    def kappa_85_40(self):
        if not self.__kappa_85_40:
            if not self.sup_kappa_sizedist or not self.sup_kappa_wavelength:
                txt = ("Make sure you define the following attributes first: \n"
                       "\t - self.sup_kappa_sizedist\n"
                       "\t - and self.sup_kappa_wavelength")
                raise ValueError(txt)

            self.__kappa_85_40, self.__growthfactor_85_40 = _hygrow.kappa_from_fofrh_and_sizedist(self.f_RH_scatt_2p_85_40, self.sup_kappa_sizedist,
                                                                                      self.sup_kappa_wavelength, 85)
        return self.__kappa_85_40

    @kappa_85_40.setter
    def kappa_85_40(self, value):
        self.__kappa_85_40 = value

    @property
    @_decorators.change_doc(_hygrow.kappa_from_fofrh_and_sizedist)
    def growth_factor_85_40(self):
        if self.__growthfactor_85_40:
            self.kappa_85_40
        return self.__growthfactor_85_40






    @property
    def sup_kappa_sizedist(self):
        return self.__sup_kappa_sizedist

    @sup_kappa_sizedist.setter
    def sup_kappa_sizedist(self, value):
        self.__kappa = None
        self.__growthfactor = None
        self.__sup_kappa_sizedist = value

    @property
    def sup_kappa_wavelength(self):
        return self.__sup_kappa_wavelength

    @sup_kappa_wavelength.setter
    def sup_kappa_wavelength(self, value):
        self.__kappa = None
        self.__growthfactor = None
        self.__sup_kappa_wavelength = value



def _concat_rules(arm_data_objs):
    """nothing here"""
    # out = arm_data_obj
    out = ArmDatasetSub(False)
    out._concat(arm_data_objs)
    return out

# def _concat_rules(arm_data_objs):
#     # create class
#     out = ArmDatasetSub(False)
#
#     # populate class with concatinated data
#     out.rh = _timeseries.concat([i.rh for i in arm_data_objs])
#     out.rh._data_periode = out._data_periode
#
#     # use time stamps from one of the variables
#     out.time_stamps = out.rh.data.index
#     return out
