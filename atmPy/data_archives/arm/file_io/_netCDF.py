import warnings as _warnings
import numpy as np
import pandas as _pd
from atmPy.general import timeseries as _timeseries
from atmPy.tools import array_tools as _arry_tools
from atmPy.tools import plt_tools as _plt_tools
import matplotlib.pylab as _plt
from atmPy.aerosols.instruments.AMS import AMS as _AMS
from atmPy.aerosols.size_distribution import sizedistribution as _sizedistribution
try:
    from netCDF4 import Dataset
except ModuleNotFoundError:
    _warnings.warn('netCDF4 not installed. You might encounter some functionality limitations.')



class Data_Quality(object):
    def __init__(self, parent, availability, availability_type, flag_info = None):
        self.parent = parent
        self.availability = availability
        self.availability_type = availability_type
        self.flag_info = flag_info

        self.__flag_matrix = None
        self.__flag_matrix_good_int_bad = None
        self.__flag_matrix_good_int_bad_either_or = None

    @property
    def flag_matrix(self):
        if type(self.__flag_matrix).__name__ == 'NoneType':
            self.__flag_matrix = self._get_flag_matrix()
        return self.__flag_matrix

    @flag_matrix.setter
    def flag_matrix(self, data):
        self.__flag_matrix = data
        self.__flag_matrix_good_int_bad = None
        self.__flag_matrix_good_int_bad_either_or = None


    @property
    def flag_matrix_good_int_bad(self):
        if type(self.__flag_matrix_good_int_bad).__name__ == 'NoneType':
            self.__flag_matrix_good_int_bad = self._get_flag_matrix_good_int_bad()
        return self.__flag_matrix_good_int_bad

    @property
    def flag_matrix_good_int_bad_either_or(self):
        if type(self.__flag_matrix_good_int_bad_either_or).__name__ == 'NoneType':
            self.__flag_matrix_good_int_bad_either_or = self._get_flag_matrix_good_int_bad_either_or()
        return self.__flag_matrix_good_int_bad_either_or

    def _get_flag_matrix_good_int_bad_either_or(self):
        flag_matrix_good_int_bad_either_or = self.flag_matrix_good_int_bad.copy()
        flag_matrix_good_int_bad_either_or['intermediate'][flag_matrix_good_int_bad_either_or['bad'] == True] = False
        return flag_matrix_good_int_bad_either_or


    def _get_flag_matrix_good_int_bad(self):
        intermediate_flags = self.flag_info[self.flag_info.quality == 'Indeterminate'].index
        bad_flags = self.flag_info[self.flag_info.quality == 'Bad'].index
        flag_matrix_good_int_bad = _pd.DataFrame(index=self.availability.index, dtype=bool)
        flag_matrix_good_int_bad['intermediate'] = _pd.Series(self.flag_matrix.loc[:, intermediate_flags].sum(axis=1).astype(bool), dtype=bool)
        flag_matrix_good_int_bad['bad'] = _pd.Series(self.flag_matrix.loc[:, bad_flags].sum(axis=1).astype(bool), dtype=bool)
        flag_matrix_good_int_bad['good'] = ~ (flag_matrix_good_int_bad['intermediate'] | flag_matrix_good_int_bad['bad'])
        return flag_matrix_good_int_bad


    def _get_flag_matrix(self):
        av = self.availability
        flag_matrix = np.zeros((av.shape[0], self.flag_info.shape[0]))
        good = np.zeros(flag_matrix.shape[0])
        for e,flag in enumerate(av.iloc[:,0]):
            b_l = list('{:0{}b}'.format(flag, self.flag_info.shape[0]))
            b_l.reverse()
            bs = np.array(b_l).astype(bool)
            flag_matrix[e] = bs
            if bs.sum() == 0:
                good[e] = 1

        flag_matrix = _pd.DataFrame(flag_matrix, dtype = bool,
                                    # index = self.parent.size_distribution.data.index,
                                    index= av.index,
                                    columns=self.flag_info.index)
        flag_matrix[0] = _pd.Series((good != 0), index = av.index)
        return flag_matrix

    def plot_stacked_bars(self, which='or', ax = None, resample=(6, 'H'), width=0.25, lw=0, show_missing = None, label = 'short', colormap = None, kwargs_leg = {}):
        """
        Args:
            which:
            ax: mpl axis instance
            resample:
            width:
            lw:
            show_missing: float
                This adds a gray background up the show_missing. Therefor show_missing should be the expacted number of points!
            kwargs_leg: dict
                kwargs passed to legend, e.g. loc, title
            label: string
                if 'short' bits are given in legend
                if 'long' flag descirption is given

        Returns:
            figur, axis
        """
        if which == 'or':
            fmrs = self.flag_matrix_good_int_bad_either_or
        elif which == 'and':
            fmrs = self.flag_matrix_good_int_bad
        elif which == 'all':
            fmrs = self.flag_matrix
        else:
            raise ValueError('{} is not an option for which'.format(which))

        if resample:
            fmrs = self._get_resampled(fmrs, resample)

        if not ax:
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()



        bars = []
        labels = []

        if which == 'all':
            if not colormap:
                colormap = _plt.cm.Accent

            bottom = np.zeros(fmrs.shape[0])
            for e,flag in enumerate(fmrs.columns):
                # b = a.bar(fmrs.index, fmrs[flag], width=width, color=_plt_tools.color_cycle[e % len(_plt_tools.color_cycle)], linewidth=lw, edgecolor='w')
                b = a.bar(fmrs.index, fmrs[flag], bottom = bottom, width=width, color=colormap(e/fmrs.columns.max()), linewidth=lw, edgecolor='w')
                bottom += fmrs[flag].values
                bars.append(b)
                if label == 'short':
                    labels.append(flag)
                elif label == 'long':
                    if flag == 0:
                        label_st = 'no flag'
                    else:
                        label_st = self.flag_info.loc[flag, 'description']
                    labels.append(label_st)
                else:
                    raise ValueError()
            # if show_missing:
            # a.legend(bars, labels, title = 'flag')
        else:
            b_g = a.bar(fmrs.index, fmrs['good'], width=width, color=_plt_tools.color_cycle[2], linewidth=lw, edgecolor='w')
            b_i = a.bar(fmrs.index, fmrs['intermediate'], width=width, bottom=fmrs['good'], color=_plt_tools.color_cycle[0], linewidth=lw, edgecolor='w')
            b_b = a.bar(fmrs.index, fmrs['bad'], width=width, bottom=fmrs['good'] + fmrs['intermediate'], color=_plt_tools.color_cycle[1], linewidth=lw, edgecolor='w')
            bars = [b_b, b_i, b_g]
            labels = ['bad', 'intermediate', 'good']
            # a.legend((b_b, b_i, b_g), ('bad', 'intermediate', 'good'))

        if show_missing:
            x = self.availability.index
            y = np.ones(x.shape) * show_missing
            cg = 0.9
            fb = a.fill_between(x, y, color=[cg, cg, cg])
            bars.append(fb)
            labels.append('missing')

        a.legend(bars, labels, **kwargs_leg)
        a.set_ylabel('Number of data points')
        a.set_xlabel('Timestamp')
        f.autofmt_xdate()
        return f, a

    def plot(self, which='or', ax = None, resample=(6, 'H'), show_missing = None, label = 'short', colormap = None, kwargs_leg = {}):
        """
        Args:
            which:
            ax: mpl axis instance
            resample:
            width:
            lw:
            show_missing: float
                This adds a gray background up the show_missing. Therefor show_missing should be the expacted number of points!
            kwargs_leg: dict
                kwargs passed to legend, e.g. loc, title
            label: string
                if 'short' bits are given in legend
                if 'long' flag descirption is given

        Returns:
            figur, axis
        """
        if which == 'or':
            fmrs = self.flag_matrix_good_int_bad_either_or
        elif which == 'and':
            fmrs = self.flag_matrix_good_int_bad
        elif which == 'all':
            fmrs = self.flag_matrix
        else:
            raise ValueError('{} is not an option for which'.format(which))

        if resample:
            fmrs = self._get_resampled(fmrs, resample)

        if not ax:
            f, a = _plt.subplots()
        else:
            a = ax
            f = a.get_figure()



        bars = []
        labels = []

        if which == 'all':
            if not colormap:
                colormap = _plt.cm.Accent

            for e,flag in enumerate(fmrs.columns):
                # b = a.bar(fmrs.index, fmrs[flag], width=width, color=_plt_tools.color_cycle[e % len(_plt_tools.color_cycle)], linewidth=lw, edgecolor='w')
                g, = a.plot(fmrs.index, fmrs[flag], color=colormap(e/fmrs.columns.max()))
                if label == 'short':
                    labels.append(flag)
                elif label == 'long':
                    if flag == 0:
                        label_st = 'no flag'
                    else:
                        label_st = self.flag_info.loc[flag, 'description']
                    labels.append(label_st)
                else:
                    raise ValueError()
            # if show_missing:
            # a.legend(bars, labels, title = 'flag')
        else:
            b_g = a.plot(fmrs.index, fmrs['good'],         color=_plt_tools.color_cycle[2])
            b_i = a.plot(fmrs.index, fmrs['intermediate'], color=_plt_tools.color_cycle[0])
            b_b = a.plot(fmrs.index, fmrs['bad'],          color=_plt_tools.color_cycle[1])
            labels = ['bad', 'intermediate', 'good']
            # a.legend((b_b, b_i, b_g), ('bad', 'intermediate', 'good'))

        if show_missing:
            x = self.availability.index
            y = np.ones(x.shape) * show_missing
            cg = 0.9
            fb = a.fill_between(x, y, color=[cg, cg, cg])
            bars.append(fb)
            labels.append('missing')

        a.legend(bars, labels, **kwargs_leg)
        a.set_ylabel('Number of data points')
        a.set_xlabel('Timestamp')
        f.autofmt_xdate()
        return f, a


    @staticmethod
    def _get_resampled(which, period=(6, 'H')):
        return which.resample(period, label='left').sum()

class ArmDataset(object):
    def __init__(self, fname, data_quality = 'good', data_quality_flag_max = None, error_bad_file = True):
        # self._data_period = None
        self._error_bad_file = error_bad_file
        self._fname = fname
        if fname:
            self.netCDF = Dataset(fname)
            self.data_quality_flag_max = data_quality_flag_max
            self.data_quality = data_quality
            if type(self.time_stamps).__name__ != 'OverflowError':
                self._parsing_error = False
                self._parse_netCDF()
            else:
                self._parsing_error = True



    def _concat(self, arm_data_objs, close_gaps = True):
        for att in self._concatable:
            first_object = getattr(arm_data_objs[0], att)
            which_type = type(first_object).__name__
            data_period = first_object._data_period
            if which_type == 'TimeSeries_2D':
                value = _timeseries.TimeSeries_2D(_pd.concat([getattr(i, att).data for i in arm_data_objs]))
            elif which_type == 'TimeSeries':
                value = _timeseries.TimeSeries(
                    _pd.concat([getattr(i, att).data for i in arm_data_objs]))
            elif which_type == 'AMS_Timeseries_lev01':
                value = _AMS.AMS_Timeseries_lev01(
                    _pd.concat([getattr(i, att).data for i in arm_data_objs]))
            elif which_type == 'SizeDist_TS':
                # value = _AMS.AMS_Timeseries_lev01(pd.concat([getattr(i, att).data for i in arm_data_objs]))
                data = _pd.concat([getattr(i, att).data for i in arm_data_objs])
                value = _sizedistribution.SizeDist_TS(data, getattr(arm_data_objs[0], att).bins, 'dNdlogDp',
                                                      ignore_data_gap_error=True,)
            elif which_type == 'TimeSeries_3D':
                value = _timeseries.TimeSeries_3D(_pd.concat([getattr(i, att).data for i in arm_data_objs]))
            else:
                raise TypeError(
                    '%s is not an allowed type here (TimeSeries_2D, TimeSeries)' % which_type)

            if hasattr(first_object, 'availability'):
                try:
                    avail_concat = _pd.concat([getattr(i, att).availability.availability for i in arm_data_objs])
                    avail = Data_Quality(None, avail_concat, None , first_object.flag_info)
                    value.availability = avail
                except:
                    _warnings.warn('availability could not be concatinated make sure you converted it to a pandas frame at some point!')
            value._data_period = data_period
            if close_gaps:
                setattr(self, att, value.close_gaps())
            else:
                setattr(self, att, value)


    @property
    def time_stamps(self):
        if '__time_stamps' in dir(self):
            return self.__time_stamps
        else:
            bt = self.netCDF.variables['base_time']
            toff = self.netCDF.variables['time_offset']
            try:
                self.__time_stamps = _pd.to_datetime(0) + _pd.to_timedelta(bt[:].flatten()[0], unit ='s') + _pd.to_timedelta(toff[:], unit ='s')
            except OverflowError as e:
                txt = str(e) + ' This probably means the netcdf file is badly shaped.'
                if self._error_bad_file:
                    raise OverflowError(txt + ' Consider setting the kwarg error_bad_file to False.')
                else:
                    _warnings.warn(txt)
                return e
            self.__time_stamps.name = 'Time'
            # self._time_offset = (60, 'm')
            if self._time_offset:
                self.__time_stamps += np.timedelta64(int(self._time_offset[0]), self._time_offset[1])
        return self.__time_stamps

    @time_stamps.setter
    def time_stamps(self,timesamps):
        self.__time_stamps = timesamps

    def _data_quality_control(self):
        return


    def _read_variable(self, variable, reverse_qc_flag = False):
        """Reads the particular variable and replaces all masked data with NaN.
        Note, if quality flag is given only values larger than the quality_control
        variable are replaced with NaN.

        Parameters
        ----------
        variable: str
            Variable name as devined in netCDF file
        reverse_qc_flag: bool or int
            Set to the number of bits, when reversion is desired
            If the indeterminate (patchy) bits are on the wrong end of the qc bit
            string it might make sense to reverte the bit string.

        Returns
        -------
        ndarray

        Examples
        --------
        self.temp = self.read_variable(ti"""
        try:
            var = self.netCDF.variables[variable]
        except:
            var = self.netCDF.variables[variable]
        data = var[:]

        variable_qc = "qc_" + variable

        availability = np.zeros(data.shape)
        availability_type = None
        if variable_qc in self.netCDF.variables.keys():
            var_qc = self.netCDF.variables["qc_" + variable]
            data_qc = var_qc[:]
            availability = data_qc.copy()
            availability_type = 'qc'

            if reverse_qc_flag:
                if type(reverse_qc_flag) != int:
                    raise TypeError('reverse_qc_flag should either be False or of type integer giving the number of bits')
                data_qc = _arry_tools.reverse_binary(data_qc, reverse_qc_flag)
            data_qc[data_qc <= self.data_quality_flag_max] = 0
            data_qc[data_qc > self.data_quality_flag_max] = 1

            # if hasattr(self, 'data_quality_max_intermediat'):
            #     print('has it!!!!')

            if data.shape != data_qc.shape:
                dt = np.zeros(data.shape)
                dt[data_qc == 1, : ] = 1
                data_qc = dt
            data = np.ma.array(data, mask = data_qc, fill_value= -9999)


        elif 'missing_data' in var.ncattrs():
            fill_value = var.missing_data
            data = np.ma.masked_where(data == fill_value, data)
        # else:
            # print('no quality flag found')
        if type(data).__name__ == 'MaskedArray':
            try:
                data.data[data.mask] = np.nan
            except ValueError:
                data = data.astype(float)
                data.data[data.mask] = np.nan
            data = data.data
        # data.availability = availability
        # data.availability_type = availability_type
        out = {}
        out['data'] = data
        out['availability'] = availability
        out['availability_type'] = availability_type
        return out

    def _read_variable2timeseries(self, variable, column_name = False, reverse_qc_flag = False):
        """
        Reads the specified variables and puts them into a timeseries.

        Parameters
        ----------
        variable: string or list of strings
            variable names
        column_name: bool or string
            this is a chance to give unites. This will also be the y-label if data
            is plotted

        Returns
        -------
        pandas.DataFrame

        """


        if type(variable).__name__ == 'str':
            variable = [variable]

        df = _pd.DataFrame(index = self.time_stamps)
        for var in variable:
            variable_out = self._read_variable(var, reverse_qc_flag = reverse_qc_flag)
            # if var == 'ratio_85by40_Bbs_R_10um_2p':
            #     import pdb
            #     pdb.set_trace()
            df[var] = _pd.Series(variable_out['data'], index = self.time_stamps)
        if column_name:
            df.columns.name = column_name
        out = _timeseries.TimeSeries(df)
        if column_name:
            out._y_label = column_name

        out._data_period = self._data_period
        out.availability = Data_Quality(self, variable_out['availability'], variable_out['availability_type'])
        # out.availability_type = variable_out['availability_type']
        return out


    def _get_variable_info(self):
        for v in self.netCDF.variables.keys():
            var = self.netCDF.variables[v]
            print(v)
            print(var.long_name)
            print(var.shape)
            print('--------')

    def _close(self):
        self.netCDF.close()
        delattr(self, 'netCDF')

    def _parse_netCDF(self):
        self._data_quality_control()
        return
