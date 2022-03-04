from atmPy.aerosols.size_distribution import sizedistribution
from atmPy.aerosols.size_distribution import diameter_binning
import pandas as pd
from atmPy.data_archives.arm._netCDF import ArmDataset
from atmPy.data_archives.arm._netCDF import Data_Quality
import numpy as np

class ArmDatasetSub(ArmDataset):
    def __init__(self,*args, **kwargs):
        self._data_period = 2700.
        self._time_offset = (- self._data_period, 's') #trying to unify position of time stamp to beginning of sampling period

        # #qc-flags
        flag_info_dict = {1: {'description':'(DMA only) Number Concentration <300 cm^-3  OR  >30,000 cm^-3','quality':'Indeterminate'},
             2: {'description':'(DMA only) Volume Concentration <0.5 um^3/cm^3  OR  >50 um^3/cm^3','quality':'Indeterminate'},
             3: {'description':'(DMA only) Unexpected volume concentration distribution slope flag at 0.12 um and 0.10 um > 1.0','quality':'Indeterminate'},
             4: {'description':'(DMA only) Poor correlation between up and down scan distributions is < 0.9','quality':'Indeterminate'},
             5: {'description':'(APS only) Volume Concentration <0.5 um^3/cm^3  OR  >50 um^3/cm^3','quality':'Indeterminate'},
             6: {'description':'(APS only) Unexpected volume concentration distribution slope flag at 0.12 um and 0.10 um > 1.0','quality':'Indeterminate'},
             7: {'description':'(Both DMA, APS) Average residual between the measured DMA and APS size distributions in the overlap of their size ranges was greater than 0.5.','quality':'Indeterminate'},
             8: {'description':'(DMA only) Number Concentration <100 cm^-3  OR  >100,000 cm^-3','quality':'Bad'},
             9: {'description':'(DMA only) Volume Concentration <0.2 um^3/cm^3  OR  >200 um^3/cm^3','quality':'Bad'},
             10:{'description':'(DMA only) Unexpected volume concentration distribution slope flag at 0.12 um and 0.10 um > 3.0','quality':'Bad'},
             11:{'description':'(DMA only) Poor correlation between up and down scan distributions is < 0.3','quality':'Bad'},
             12:{'description':'(APS only) Volume Concentration <0.2 um^3/cm^3  OR  >200 um^3/cm^3','quality':'Bad'},
             13:{'description':'(APS only) Unexpected volume concentration distribution slope flag at 0.12 um and 0.10 um > 3.0','quality':'Bad'},
             14:{'description':'(Both DMA, APS) Average residual between the measured DMA and APS size distributions in the overlap of their size ranges was greater than 0.8.','quality':'Bad'},
             }

        flag_info = pd.DataFrame(flag_info_dict).transpose().copy()
        flag_info.index.name = 'Bit'
        self.flag_info  = flag_info

        super(ArmDatasetSub,self).__init__(*args, **kwargs)
        self._concatable = ['size_distribution']



    def _data_quality_control(self):
        self.data_quality_max_intermediat = 127
        if self.data_quality_flag_max == None:
            if self.data_quality == 'good':
                self.data_quality_flag_max = 0
            elif self.data_quality == 'patchy':
                self.data_quality_flag_max = self.data_quality_max_intermediat
            elif self.data_quality == 'bad':
                self.data_quality_flag_max = np.inf # 100000
            else:
                txt = '%s is not an excepted values for data_quality ("good", "patchy", "bad")'%(self.data_quality)
                raise ValueError(txt)

    def _parse_netCDF(self):
        super(ArmDatasetSub,self)._parse_netCDF()

        data = self._read_variable('number_concentration_DMA_APS')
        df = pd.DataFrame(data['data'],
                          index = self.time_stamps)

        d = self._read_variable('diameter')['data']
        bins, colnames = diameter_binning.bincenters2binsANDnames(d[:]*1000)

        self.size_distribution = sizedistribution.SizeDist_TS(df,bins,'dNdlogDp', ignore_data_gap_error = True,
                                                              # fill_data_gaps_with = np.nan
                                                              )
        self.size_distribution._data_period = self._data_period
        self.size_distribution.flag_info = self.flag_info
        availability = pd.DataFrame(data['availability'], index = self.time_stamps)
        self.size_distribution.availability = Data_Quality(self, availability, data['availability_type'], self.flag_info)
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
