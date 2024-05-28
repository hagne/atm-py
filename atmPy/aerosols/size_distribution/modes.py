from . import sizedistribution
import matplotlib.pylab as plt
import numpy as np
import pandas as pd
from scipy import signal
from scipy import optimize
from ...tools import math_functions


def fit_normal_dist(sd, log=True, p0=[10, 180, 0.2], boundary_width = [0.1, 0.6]):
    """Fits a normal distribution to a """
    if not log:
        txt = 'sorry, this is not working right now ... programming requried'
        raise ValueError(txt)

    def find_peak_arg(x, y, start_w=0.2, tol=0.3):
        """
        Parameters
        ----------
        x: nd_array
            log10 of the diameters
        y: nd_array
            intensities (number, volume, surface ...)
        start_w: float
            some reasonalble width for a log normal distribution (in log normal)
        tol: float
            Tolerance ratio for start_w"""

        med = np.median(x[1:] - x[:-1])

        low = np.floor((start_w * (1 - tol)) / med)
        top = np.ceil((start_w * (1 + tol)) / med)
        
        if low < 1:
            low = 1
        
        widths = np.arange(low, top)
        peakind = signal.find_peaks_cwt(y, widths)
        return peakind

    def multi_gauss(x, *params, verbose=False):
        y = np.zeros_like(x)
        for i in range(0, len(params), 3):
            if verbose:
                print(len(params), i)
            amp = params[i]
            pos = params[i + 1]
            sig = params[i + 2]
            y = y + math_functions.gauss(x, amp, pos, sig)
        return y

    out = {}

    x = sd.data.columns.values
    y = sd.data.iloc[0, :].values
    x = x[~ np.isnan(y)]
    y = y[~ np.isnan(y)]

    if len(x) == 0:
        return False

    if log:
        x_orig = x.copy()  # we have to keep it to avoide rounding errors when doing a back and forth calculation
        x = np.log10(x)

    out['as_fitted'] = (x, y)

    start_width = 0.2
    tol = 0.9
    # width_ll = start_width * (1 - tol)
    # width_ul = start_width * (1 + tol)
    peak_args = find_peak_arg(x, y, start_w=start_width, tol=tol)
    out['peak_args'] = peak_args

    param = []
    bound_l = []
    bound_h = []
    for pa in peak_args:
        # amp
        #         print('amp: ', y[pa])
        param.append(y[pa])
        #         bound_l.append(-np.inf)
        bound_l.append(0)
        bound_h.append(np.inf)

        # pos
        #         print('pos: ', 10**x[pa])
        param.append(x[pa])
        bound_l.append(x[pa] - 0.1)
        bound_h.append(x[pa] + 0.1)
        #         bound_l.append(0)
        #         bound_h.append(x[-1])

        # sig
        param.append(start_width)
        bound_l.append(boundary_width[0])
        bound_h.append(boundary_width[1])

    try:
        param, cm = optimize.curve_fit(multi_gauss, x, y, p0=param, bounds=(bound_l, bound_h))
    except RuntimeError:
        print('Runtime error!')
        return False
    y_fit = multi_gauss(x, *param)

    param = param.reshape(len(peak_args), 3).transpose()
    param_df = pd.DataFrame(param.transpose(), columns=['amp', 'pos', 'sig'])
    out['fit_res_param_pre'] = param_df.copy()

    #### mode attribution
    gaus = pd.DataFrame(index=x)
    for idx in param_df.index:
        gaus[idx] = pd.Series(math_functions.gauss(x, param_df.loc[idx, 'amp'], param_df.loc[idx, 'pos'], param_df.loc[idx, 'sig']), index=x)

    sum_peaks = gaus.sum(axis=1)

    gaus_rel = gaus.copy()
    for col in gaus.columns:
        gaus_rel[col] = gaus.loc[:, col] / sum_peaks

    dist_by_type = gaus_rel.copy()
    for col in gaus_rel.columns:
        dist_by_type[col] = gaus_rel.loc[:, col] * y

    #### fix x axis back to diameter
    if log:
        param[1] = 10 ** param[1]
        #         x = 10 ** x
        dist_by_type.index = x_orig

    param_df.index.name = 'peak'
    out['fit_res_param'] = param_df
    out['fit_res'] = pd.DataFrame(y_fit, index=x)
    dist_by_type.index.name = 'bin_center(nm)'
    dist_by_type.columns.name = 'peak'
    out['dist_by_type'] = dist_by_type
    return out


class FitRes(object):
    def __init__(self, fitres):
        self.fitres = fitres

    def plot(self, ax = None, scale = 2000, **kwargs):
        post = self.fitres['pos'].copy()
        post = np.log10(post)
        p_min = post.min()
        p_max = post.max() - p_min
        # cols = list(((post - p_min) / p_max).values)
        # cols = plt.cm.Accent(cols)

        if not ax:
            f ,a = plt.subplots()
        else:
            a = ax
            f = a.get_figure()
        a.scatter(self.fitres['pos'].index, self.fitres['pos'], s = self.fitres['area_rel'] * scale,
                  # color = cols,
                  **kwargs)
        # g = a.get_lines()[-1]
        # g.set_markersize(fit_res_all['area_rel'].values)
        a.set_yscale('log')
        a.set_xlim(self.fitres['pos'].index[0], self.fitres['pos'].index[-1])
        f.autofmt_xdate()
        return a


class ModeAnalysis(object):
    def __init__(self, sizedist):
        self._parent = sizedist

        self.__size_dist_aiken = None
        self.__size_dist_accu = None
        self.__size_dist_coarse = None
        self.__mode_fit_results = None

    @property
    def fit_results(self):
        if not self.__mode_fit_results:
            self.find_modes()
        return self.__mode_fit_results

    @property
    def size_dist_aiken(self):
        if not self.__size_dist_aiken:
            self.find_modes()
        return self.__size_dist_aiken

    @property
    def size_dist_accu(self):
        if not self.__size_dist_accu:
            self.find_modes()
        return self.__size_dist_accu

    @property
    def size_dist_coarse(self):
        if not self.__size_dist_coarse:
            self.find_modes()
        return self.__size_dist_coarse

    def find_modes(self):
        """This function will try to find different aerosol modes in sizedist.

        Parameters
        ----------
        sizedist: sizedistribution instances (SizeDist, SizeDist_TS, ...)
        """
        sizedist = self._parent
        boundary_accu_coars = 1000
        boundary_aiken_accu = 50

        sdts_aiken = sizedist.copy()
        sdts_aiken._update()
        sdts_aiken.data[:] = np.nan
        sdts_accu = sdts_aiken.copy()
        sdts_accu._update()
        sdts_coarse = sdts_aiken.copy()
        sdts_coarse._update()

        fit_res_all = pd.DataFrame(columns=['amp', 'pos', 'sig', 'area'])

        for which in sizedist.data.index:
            #     which = sizedist.data.index[21]
            sd = sizedistribution.SizeDist(pd.DataFrame(sizedist.data.loc[which ,:]).transpose(), sizedist.bins, sizedist.distributionType)
            #     sd = sd.convert2dVdlogDp()
            out_f = fit_normal_dist(sd)
            if not out_f:
                continue
            # res_dict[which] = out_f

            dist_by_type = out_f['dist_by_type'].copy()
            fit_res_param = out_f['fit_res_param']

            coarse = fit_res_param.loc[: ,'pos'] > boundary_accu_coars
            aiken = fit_res_param.loc[: ,'pos'] < boundary_aiken_accu
            accu = np.logical_and(fit_res_param.loc[: ,'pos'] <= boundary_accu_coars, fit_res_param.loc[: ,'pos'] >= boundary_aiken_accu)

            df = pd.DataFrame(index = dist_by_type.index)
            df['aiken'] = dist_by_type.iloc[: ,aiken.values].sum(axis = 1)
            df['accu'] = dist_by_type.iloc[: ,accu.values].sum(axis = 1)
            df['coarse'] = dist_by_type.iloc[: ,coarse.values].sum(axis = 1)
            sdts_aiken.data.loc[which ,:] = df.loc[: ,'aiken']
            sdts_accu.data.loc[which ,:] = df.loc[: ,'accu']
            sdts_coarse.data.loc[which ,:] = df.loc[: ,'coarse']

            fit_res_param['area'] = fit_res_param['amp'] * fit_res_param['sig'] * np.sqrt(2 * np.pi)
            fit_res_param['area_rel'] = fit_res_param['area'] / fit_res_param['area'].sum()
            fit_res_param.index = [which] * fit_res_param.shape[0]
            # fit_res_all = fit_res_all.append(fit_res_param) # deprecated
            fit_res_all = pd.concat([fit_res_all, fit_res_param], 
                                    # ignore_index=True,
                                   )
            if 0:
                f ,a = sd.plot()
                a.set_title(which)
                sdts_aiken.data.loc[which ,:].plot(ax = a)
                sdts_accu.data.loc[which ,:].plot(ax = a)
                sdts_coarse.data.loc[which ,:].plot(ax = a)
        # (sdts_coarse.data.loc[which,:] + sdts_accu.data.loc[which,:]).plot(ax = a, color = 'magenta')
        #     break
        # df

        # volumes = sdts_aiken.particle_volume_concentration.data.copy()
        # volumes.columns = ['aiken']
        # volumes['acccu'] = sdts_accu.particle_volume_concentration.data
        # volumes['coarse'] = sdts_coarse.particle_volume_concentration.data
        # volume_ratios = volumes.truediv(volumes.sum(axis = 1), axis=0)

        self.__size_dist_aiken = sdts_aiken
        self.__size_dist_accu = sdts_accu
        self.__size_dist_coarse = sdts_coarse

        self.__mode_fit_results = FitRes(fit_res_all)
