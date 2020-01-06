import pandas as _pd
import matplotlib.pylab as _plt
import numpy as _np
from matplotlib.dates import MonthLocator  as _MonthLocator
from matplotlib.dates import DateFormatter as _DateFormatter
from matplotlib.ticker import FuncFormatter as _FuncFormatter
from matplotlib.ticker import MultipleLocator as _MultipleLocator
import plt_tools as _plt_tools
from atmPy.tools import array_tools as _array_tools

class Statistics(object):
    def __init__(self, parent_ts):
        self._parent_ts = parent_ts
        self.seasonality = Seasonality(self)
        self.diurnality = Diurnality(self)

    def define_custom(self,reference = 'index', frequency = 'M', bins = None):
        self.custom = Climatology(self, reference = reference, frequency = frequency, bins = bins)


class Climatology(object):
    def __init__(self, parent_stats, reference = 'index', frequency = 'M', bins = None):
        """

        Parameters
        ----------
        parent_stats
        reference: str ['index']
            if one wants to plot agains a particular column in the timeseries rather then the index.
        frequency: str (['M'], 'H')
            originally this was written for things like seasonality, diurnality, with otherwords with anyting related
            to the datetime.
        bins: array like
            This is only used if reference column (e.g. index) is not a datetime object, e.g. if other climatologic
            indeces that are not datetime related are used.
        """
        self._bins = bins
        self._reference = reference
        self._parent_stats = parent_stats
        self._parent_ts = parent_stats._parent_ts
        self._frequency = frequency
        self._reset()
        self._timezone = None

    def _reset(self):
        self._percentiles = None


    @property
    def frequency(self):
        return self._frequency

    @property
    def percentiles(self):
        if type(self._percentiles) == type(None):
            # rs = self._parent_ts.data.resample(self.frequency, label='left',
            #                                    # convention = 'start',
            #                                    # closed = 'left'
            #                                    )
            data = self._parent_ts.data.copy()

            # if another column instead of the index is to be used as a reference
            if not self._reference == 'index':
                data.index = data[self._reference]
                data.drop(labels = [self._reference], axis = 1, inplace = True)

            if type(data.index).__name__ == 'DatetimeIndex':
                if self._timezone:
                    data.index += _pd.Timedelta(self._timezone, 'h')
                if self.frequency == 'H':
                    data.index = data.index.hour
                elif self.frequency == 'M':
                    data.index = data.index.month
            else:
                data.index = _pd.Series(data.index).apply(lambda x: self._bins[abs(self._bins - x).argmin()])

            data.sort_index(inplace=True)
            rs = data.groupby(data.index)
            def percentile(x, q):
                # print(x.columns)

                # print('-----')
                x = x.dropna()
                if x.shape[0] == 0:
                    pct = _np.nan
                elif x.shape[1] == 1:
                    pct = _np.percentile(x, q)
                elif x.shape[1] == 2:
                    weights = x.weights
                    values = x.drop(['weights'], axis = 1).iloc[:,0]
                    # print(values.)
                    # print(x.columns)
                    # print([q/100.])
                    # print('--------')
                    pct = _array_tools.weighted_quantile(values, [q/100.], sample_weight = weights)[0]
                return pct

            def average(x):
                x = x.dropna()
                weights = x.weights
                values = x.drop(['weights'], axis=1).iloc[:, 0]
                avg = _np.average(values,weights = weights)
                return avg

            def median(x):
                x = x.dropna()
                cols = list(x.columns)
                # print(cols)
                cols.pop(cols.index('weights'))
                x.sort_values(cols[0], inplace=True)
                cumsum = x.weights.cumsum()
                cutoff = x.weights.sum() / 2.0
                median = x[cols[0]][cumsum >= cutoff].iloc[0]
                return median

            def number_of_valid_data_points(x):
                x = x.dropna()
                return x.shape[0]

            out = _pd.DataFrame()
            if data.shape[1] == 1:
                out['mean'] = rs.mean().iloc[:, 0]
                out['median'] = rs.median().iloc[:,0]
            elif data.shape[1] == 2:
                if 'weights' not in data.columns:
                    raise KeyError('If two columns are given one of them must have the label "weights"')
                out['mean'] = rs.apply(average)
                out['median'] = rs.apply(median)

            perc_list = [5, 25, 50, 75, 95]
            for perc in perc_list:
                outt = rs.apply(lambda x: percentile(x, perc))
                # print(outt.shape)
                out[perc] = outt#.iloc[:, 0]
            out['n_valid'] = rs.apply(number_of_valid_data_points)
            # out.index += _np.timedelta64(1, 'D')
            self._percentiles = out
        return self._percentiles

    def plot_percentiles(self, ax=None, box_width=0.2, wisker_size=20, mean_size=10, median_size = 10 , line_width=1.5,
                         xoffset=0,
                         color=0, tickbase = 1):
        """
        This will plot the percentiles of the data ([5, 25, 50, 75, 95]. In addition, the mean (o) and median (_) are shown.

        Parameters
        ----------
        ax
        box_width
        wisker_size
        mean_size
        median_size
        line_width
        xoffset
        color
        tickbase

        Returns
        -------
        f, a, boxes, vlines, wisker_tips, mean
        """
        if type(color) == int:
            color = _plt.rcParams['axes.prop_cycle'].by_key()['color'][color]
            col = _plt_tools.colors.Color(color, model='hex')
        elif type(color) == str:
            col = _plt_tools.colors.Color(color, model='hex')
        else:
            col = _plt_tools.colors.Color(color, model='rgb')

        col.saturation = 0.3
        color_bright = col.rgb

        if ax:
            a = ax
            f = a.get_figure()
        else:
            f, a = _plt.subplots()

        boxes = []
        vlines = []
        xordinal = []
        for row in self.percentiles.iterrows():
            #             width = 10
            x = row[0] + xoffset
            xordinal.append(x)

            # box
            # y = (row[1][75] + row[1][25]) / 2
            y = row[1][25]
            height = row[1][75] - row[1][25]
            box = _plt.Rectangle((x - box_width / 2, y), box_width, height,
                                 #                         ha = 'center'
                                 )
            box.set_facecolor([1, 1, 1, 1])
            a.add_patch(box)
            boxes.append(box)
            # wiskers
            y = (row[1][95] + row[1][5]) / 2
            vl = a.vlines(x, row[1][5], row[1][95])
            vlines.append(vl)

        for b in boxes:
            b.set_linewidth(line_width)
            b.set_facecolor(color_bright)
            b.set_edgecolor(color)
            b.set_zorder(2)

        for vl in vlines:
            vl.set_color(color)
            vl.set_linewidth(line_width)
            vl.set_zorder(1)

        # wm_lw = 2
        wisker_tips = []
        if wisker_size:
            g, = a.plot(xordinal, self.percentiles[5], ls='')
            wisker_tips.append(g)

            g, = a.plot(xordinal, self.percentiles[95], ls='')
            wisker_tips.append(g)

        for wt in wisker_tips:
            wt.set_markeredgewidth(line_width)
            wt.set_color(color)
            wt.set_markersize(wisker_size)
            wt.set_marker('_')

        mean = None
        if mean_size:
            g, = a.plot(xordinal, self.percentiles['mean'], ls='')
            g.set_marker('o')
            g.set_markersize(mean_size)
            g.set_zorder(20)
            g.set_markerfacecolor('None')
            g.set_markeredgewidth(line_width)
            g.set_markeredgecolor(color)
            mean = g

        median = None
        if median_size:
            g, = a.plot(xordinal, self.percentiles['median'], ls='')
            g.set_marker('_')
            g.set_markersize(median_size)
            g.set_zorder(20)
            g.set_markeredgewidth(line_width)
            g.set_markeredgecolor(color)
            median = g

        # a.xaxis.set_major_locator(_MonthLocator())
        # a.xaxis.set_major_formatter(_DateFormatter('%b'))
        try:
            a.set_ylim(_np.nanmin(self.percentiles.drop(['n_valid'], axis=1)),
                       _np.nanmax(self.percentiles.drop(['n_valid'], axis=1)))
        except:
            pass

        try:
            a.set_xlim(self.percentiles.index.min(), self.percentiles.index.max())
        except:
            pass

        mjl = _MultipleLocator(tickbase)
        a.xaxis.set_major_locator(mjl)




        # a.relim()
        # a.autoscale_view(tight=True)
        # f.autofmt_xdate()
        return f, a, boxes, vlines, wisker_tips, mean, median

            # def plot_percentiles(self, ax=None, box_width=10, wisker_size=20, mean_size = 10, line_width = 1.5, xoffset = 0, color=0):
    #     """
    #
    #     Parameters
    #     ----------
    #     ax
    #     box_width
    #     wisker_size
    #     mean_size
    #     color
    #
    #     Returns
    #     -------
    #     f, a, boxes, vlines, wisker_tips, mean
    #     """
    #     if type(color) == int:
    #         color = _plt.rcParams['axes.prop_cycle'].by_key()['color'][color]
    #         col = _plt_tools.colors.Color(color, model='hex')
    #     else:
    #         col = _plt_tools.colors.Color(color, model='rgb')
    #
    #     col.saturation = 0.3
    #     color_bright = col.rgb
    #
    #     if ax:
    #         a = ax
    #         f = a.get_figure()
    #     else:
    #         f, a = _plt.subplots()
    #
    #     boxes = []
    #     vlines = []
    #     xordinal = []
    #     for row in self.percentiles.iterrows():
    #         #             width = 10
    #         x = row[0].toordinal() + xoffset
    #         xordinal.append(x)
    #
    #         # box
    #         y = (row[1][75] + row[1][25]) / 2
    #         height = row[1][75] - row[1][25]
    #         box = _plt.Rectangle((x - box_width / 2, y), box_width, height,
    #                              #                         ha = 'center'
    #                              )
    #         box.set_facecolor([1, 1, 1, 1])
    #         a.add_patch(box)
    #         boxes.append(box)
    #         # wiskers
    #         y = (row[1][95] + row[1][5]) / 2
    #         vl = a.vlines(x, row[1][5], row[1][95])
    #         vlines.append(vl)
    #
    #     for b in boxes:
    #         b.set_linewidth(line_width)
    #         b.set_facecolor(color_bright)
    #         b.set_edgecolor(color)
    #         b.set_zorder(2)
    #
    #     for vl in vlines:
    #         vl.set_color(color)
    #         vl.set_linewidth(line_width)
    #         vl.set_zorder(1)
    #
    #     # wm_lw = 2
    #     wisker_tips = []
    #     if wisker_size:
    #         g, = a.plot(xordinal, self.percentiles[5], ls='')
    #         wisker_tips.append(g)
    #
    #         g, = a.plot(xordinal, self.percentiles[95], ls='')
    #         wisker_tips.append(g)
    #
    #     for wt in wisker_tips:
    #         wt.set_markeredgewidth(line_width)
    #         wt.set_color(color)
    #         wt.set_markersize(wisker_size)
    #         wt.set_marker('_')
    #
    #     mean = None
    #     if mean_size:
    #         g, = a.plot(xordinal, self.percentiles['mean'], ls = '')
    #         g.set_marker('o')
    #         g.set_markersize(mean_size)
    #         g.set_zorder(20)
    #         g.set_markerfacecolor('None')
    #         g.set_markeredgewidth(line_width)
    #         g.set_markeredgecolor(color)
    #         mean = g
    #
    #     a.xaxis.set_major_locator(_MonthLocator())
    #     a.xaxis.set_major_formatter(_DateFormatter('%b'))
    #
    #     a.set_ylim(_np.nanmin(self.percentiles.drop(['n_valid'], axis=1)), _np.nanmax(self.percentiles.drop(['n_valid'], axis=1)))
    #     a.set_xlim(self.percentiles.index.min().toordinal(), self.percentiles.index.max().toordinal())
    #     # a.relim()
    #     # a.autoscale_view(tight=True)
    #     # f.autofmt_xdate()
    #     return f, a, boxes, vlines, wisker_tips, mean

class Seasonality(Climatology):
    def plot_percentiles(self, *args, **kwargs):
        out = super().plot_percentiles(*args, **kwargs)
        a = out[1]
        def num2month(pos, num):
            month = ['', 'J', 'F', 'M', 'A', 'M', 'J', 'J', 'A', 'S', 'O', 'N', 'D', '', '']
            return month[num]

        if self.frequency == 'M':
            a.xaxis.set_major_formatter(_FuncFormatter(num2month))
            a.set_xlim(0.5, 12.5)
        a.set_xlabel('Month of year')
        return out

class Diurnality(Climatology):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, frequency='H')

    def plot_percentiles(self, *args, **kwargs):
        out = super().plot_percentiles(*args, **kwargs)
        a = out[1]

        if self.frequency == 'H':
            a.set_xlim(-0.5,23.5)

        a.set_xlabel('Hour of day')
        return out

    @property
    def timezone(self):
        return self._timezone

    @timezone.setter
    def timezone(self, value):
        self._reset()
        self._timezone = value


