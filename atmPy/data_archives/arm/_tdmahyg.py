import numpy as np
import pandas as pd

from atmPy.aerosols.physics import hygroscopic_growth as hg
from atmPy.general import timeseries
from atmPy.data_archives.arm._netCDF import ArmDataset


class ArmDatasetSub(ArmDataset):
    def __init__(self,*args, **kwargs):
        self._data_period = 2700.
        self._time_offset = (- self._data_period, 's')
        super(ArmDatasetSub,self).__init__(*args, **kwargs)
        self._concatable = ['RH_interDMA', 'hyg_distributions']
        self.__kappa_values = None


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
        "returns a dictionary, with panels in it"
        super(ArmDatasetSub,self)._parse_netCDF()


        size_bins = self._read_variable('size_bins') * 1000
        df = pd.DataFrame(self._read_variable('RH_interDMA'), index = self.time_stamps, columns=size_bins)
        df.columns.name = 'size_bin_center_nm'
        self.RH_interDMA = timeseries.TimeSeries(df)
        self.RH_interDMA._data_period = self._data_period

        data = self._read_variable('hyg_distributions')
        growthfactors = self._read_variable('growthfactors')
        data = pd.Panel(data, items= self.time_stamps, major_axis = size_bins, minor_axis = growthfactors)
        data.major_axis.name = 'size_bin_center_nm'
        data.minor_axis.name = 'growthfactors'
        self.hyg_distributions = timeseries.TimeSeries_3D(data)
        self.hyg_distributions._data_period = self._data_period

    def plot_all(self):
        self.hyg_distributions.plot(yaxis=2, sub_set=5)
        self.RH_interDMA.plot()

    @property
    def mean_growth_factor(self):
        """Calculates the mean growthfactor of the particular size bin."""
        if '__mean_growth_factor' not in dir(self):


            def mean_linewise(gf_dist):
                growthfactors = self.hyg_distributions.data.minor_axis.values
                # meanl = ((gf_dist[~ gf_dist.mask] * np.log10(growthfactors[~ gf_dist.mask])).sum()/gf_dist[~gf_dist.mask].sum())
                meanl = ((gf_dist[~ np.isnan(gf_dist)] * np.log10(growthfactors[~ np.isnan(gf_dist)])).sum()/gf_dist[~np.isnan(gf_dist)].sum())
                stdl = np.sqrt((gf_dist[~ np.isnan(gf_dist)] * (np.log10(growthfactors[~ np.isnan(gf_dist)]) - meanl)**2).sum()/gf_dist[~np.isnan(gf_dist)].sum())
                return np.array([10**meanl,stdl])
            data = self.hyg_distributions.data
            allmeans = timeseries.TimeSeries_3D(pd.Panel(items=data.items, major_axis=data.major_axis, minor_axis= ['mean', 'std_log']))
            for i,time in enumerate(data.values):
                for e,size in enumerate(time):
                    allmeans.data.iloc[i,e] = mean_linewise(size)
            self.__mean_growth_factor = allmeans
            self.__mean_growth_factor._data_period = self._data_period
        return self.__mean_growth_factor

    @property
    def kappa_values(self):
        if not self.__kappa_values:
            # RH =
            kappa_values = hg.kappa_simple(self.mean_growth_factor.data.values[:,:,0],self.RH_interDMA.data.values, inverse = True)
            kappa_values = pd.DataFrame(kappa_values,columns=self.mean_growth_factor.data.major_axis, index = self.mean_growth_factor.data.items)
            self.__kappa_values = timeseries.TimeSeries_2D(kappa_values)
            # self.plottable.append('kappa_values')
            self.__kappa_values._data_period = self._data_period
        return self.__kappa_values

    @kappa_values.setter
    def kappa_values(self, value):
        self.__kappa_values = value

def _concat_rules(arm_data_objs):
    """nothing here"""
    # out = arm_data_obj
    out = ArmDatasetSub(False)
    out._concat(arm_data_objs)
    return out

# def _concat_rules(arm_data_objs):
#     out = ArmDatasetSub(False)
#     out.RH_interDMA = timeseries.TimeSeries(pd.concat([i.RH_interDMA.data for i in arm_data_objs]))
#     out.RH_interDMA._data_period = out._data_period
#
#     out.hyg_distributions = timeseries.TimeSeries_3D(pd.concat([i.hyg_distributions.data for i in arm_data_objs]))
#     out.hyg_distributions._data_period = out._data_period
#
#     out.time_stamps = out.RH_interDMA.data.index
#     return out
