from atmPy.general import timeseries as _timeseries
from atmPy.data_archives.arm import _netCDF


class ArmDatasetSub(_netCDF.ArmDataset):
    def __init__(self,*args, **kwargs):
        self._data_period = None
        self._time_offset = (- self._data_period, 's')
        super(ArmDatasetSub,self).__init__(*args, **kwargs)

        ####
        # for properties
        self.__mean_growth_factor  = None
        self._concatable = ['rh']


    def _data_quality_control(self):
        #######
        ## Define what is good, patchy or bad data
        ## delete if not quality flags exist
        ## settings are ignored if no quality flags are profided (e.g. noaaaos)

        if self.data_quality_flag_max == None:
            if self.data_quality == 'good':
                self.data_quality_flag_max = 0
            elif self.data_quality == 'patchy':
                self.data_quality_flag_max = 127
            elif self.data_quality == 'bad':
                self.data_quality_flag_max = 100000
            else:
                txt = '%s is not an excepted values for data_quality ("good", "patchy", "bad")'%(self.data_quality)
                raise ValueError(txt)

    def _parse_netCDF(self):
        super(ArmDatasetSub,self)._parse_netCDF()
        self.rh = self._read_variable2timeseries(['rh_60m', 'rh_60m'], column_name='Relative Humidity (%)')

    def plot_all(self):
        self.rh.plot()

    @property
    def mean_growth_factor(self):
        """do something more with the data"""
        if not self.__mean_growth_factor:
            self.__mean_growth_factor = 'bla'
        return self.__mean_growth_factor


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
