from atmPy.general import timeseries as _timeseries
from atmPy.data_archives.arm import _netCDF


class ArmDatasetSub(_netCDF.ArmDataset):
    def __init__(self,*args, **kwargs):
        self._data_period = 60.
        self._time_offset = (- self._data_period, 's')
        super(ArmDatasetSub,self).__init__(*args, **kwargs)
        ## Define what is good, patchy or bad data

        # self._parse_netCDF()

    def _parse_netCDF(self):
        super(ArmDatasetSub,self)._parse_netCDF()
        # self._data_quality_control()
        self.relative_humidity = self._read_variable2timeseries(['rh_25m', 'rh_60m'], column_name='Relative Humidity (%)')
        self.temperature = self._read_variable2timeseries(['temp_25m', 'temp_60m'], column_name='Temperature ($^{\circ}$C)')
        self.vapor_pressure = self._read_variable2timeseries(['vap_pres_25m', 'vap_pres_60m'], column_name='Vapor pressure (kPa)')

    def _data_quality_control(self):
        if self.data_quality_flag_max == None:
            if self.data_quality == 'good':
                self.data_quality_flag_max = 0
            elif self.data_quality == 'patchy':
                self.data_quality_flag_max = 0
            elif self.data_quality == 'bad':
                self.data_quality_flag_max = 100000
            else:
                txt = '%s is not an excepted values for data_quality ("good", "patchy", "bad")'%(self.data_quality)
                raise ValueError(txt)

    def plot_all(self):
        self.relative_humidity.plot()
        self.temperature.plot()
        self.vapor_pressure.plot()


def _concat_rules(arm_data_objs):
    # create class
    out = ArmDatasetSub(False)

    # populate class with concatinated data
    out.relative_humidity = _timeseries.concat([i.relative_humidity for i in arm_data_objs])
    out.relative_humidity._data_period = out._data_period
    out.temperature = _timeseries.concat([i.temperature for i in arm_data_objs])
    out.temperature._data_period = out._data_period
    out.vapor_pressure = _timeseries.concat([i.vapor_pressure for i in arm_data_objs])
    out.vapor_pressure._data_period = out._data_period

    # use time stamps from one of the variables
    out.time_stamps = out.relative_humidity.data.index
    return out
