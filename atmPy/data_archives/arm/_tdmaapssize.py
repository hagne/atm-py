from atmPy.aerosols.size_distribution import sizedistribution
from atmPy.aerosols.size_distribution import diameter_binning
import pandas as pd
from atmPy.data_archives.arm._netCDF import ArmDataset


class ArmDatasetSub(ArmDataset):
    def __init__(self,*args, **kwargs):
        self._data_period = 2700.
        self._time_offset = (- self._data_period, 's')
        super(ArmDatasetSub,self).__init__(*args, **kwargs)
        self._concatable = ['size_distribution']


    def _data_quality_control(self):
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

        df = pd.DataFrame(self._read_variable('number_concentration_DMA_APS'),
                          index = self.time_stamps)

        d = self._read_variable('diameter')
        bins, colnames = diameter_binning.bincenters2binsANDnames(d[:]*1000)

        self.size_distribution = sizedistribution.SizeDist_TS(df,bins,'dNdlogDp')
        self.size_distribution._data_period = self._data_period
        # conc = self.read_variable()

    def plot_all(self):
        self.size_distribution.plot()

def _concat_rules(arm_data_objs):
    """nothing here"""
    # out = arm_data_obj
    out = ArmDatasetSub(False)
    out._concat(arm_data_objs)
    return out

# def _concat_rules(files):
#     out = ArmDatasetSub(False)
#     data = pd.concat([i.size_distribution.data for i in files])
#     out.size_distribution = sizedistribution.SizeDist_TS(data,files[0].size_distribution.bins,'dNdlogDp')
#     out.size_distribution._data_period = out._data_period
#     return out
