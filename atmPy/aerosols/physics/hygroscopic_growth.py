import numpy as _np
from scipy.optimize import fsolve as _fsolve
import pandas as _pd
import atmPy.general.timeseries as _timeseries
import warnings as _warnings

def kappa_simple(k, RH, refractive_index = None, inverse = False):
    """Returns the growth factor as a function of kappa and RH.
    This function is based on the simplified model introduced by Rissler et al. (2006).
    It ignores the curvature effect and is therefore independend of the exact particle diameter.
    Due to this simplification this function is valid only for particles larger than 100 nm
    Patters and Kreidenweis (2007).

    Petters, M. D., & Kreidenweis, S. M. (2007). A single parameter representation of hygroscopic
    growth and cloud condensation nucleus activity, 1961-1971. doi:10.5194/acp-7-1961-2007

    Rissler, J., Vestin, A., Swietlicki, E., Fisch, G., Zhou, J., Artaxo, P., & Andreae, M. O. (2005).
    Size distribution and hygroscopic properties of aerosol particles from dry-season biomass burning
    in Amazonia. Atmospheric Chemistry and Physics Discussions, 5(5), 8149-8207. doi:10.5194/acpd-5-8149-2005

    latex expression: $\sqrt[3]{1 + \kappa \cdot \frac{RH}{100 - RH}}$

    Arguments
    ---------
    k: float
        kappa value between 0 (no hygroscopicity) and 1.4 (NaCl; most hygroscopic
        material found in atmospheric
        aerosols)
        If inverse is True than this is the growth factor.
    RH: float
        Relative humidity -> between 0 and 100 (i don't think this simplified
        version will be correct at higher RH?!?)
    inverse: bool.
        if function is inversed, which means instead of the growth factor is
        calculated as a function of the kappa value, it is the other way around.
    Returns
    -------
    float: The growth factor of the particle.
    if refractive_index is given a further float is returned which gives the new refractiv index of the grown particle
    """

    if not inverse:
        if _np.any(k > 1.4) or _np.any(k < 0):
            # txt = '''The kappa value has to be between 0 and 1.4.'''
            # txt = txt + '\nk[k>1.4]:\n%s'%(k[k>1.4]) +'\nk[k<0]:\n%s'%(k[k<0])
            # print(txt)
            # raise ValueError(txt)
            if _np.any(k > 1.4):
                _warnings.warn('There are kappa values lareger than 1.4 in this dataset!')
            if _np.any(k < 0):
                _warnings.warn('There are kappa values smaller than 0 in this dataset!')

    if _np.any(RH > 100) or _np.any(RH < 0):
        txt = """RH has to be between 0 and 100"""
        raise ValueError(txt)

    if inverse:
        kappa = lambda gf,RH: (gf**3 - 1) / (RH / (100 - RH))
        out = kappa(k,RH)
    else:
        growths_factor = lambda k,RH: (1 + (k * (RH/(100 - RH))))**(1/3.)
        out = growths_factor(k,RH)

    # adjust index of refraction
    if not inverse:
        if _np.any(refractive_index):
            nw = 1.33
            # n_mix = lambda n,gf: (n + ((gf-1)*nw))/gf
            n_mix = lambda n, gf: (n + (nw * (gf ** 3 - 1))) / gf ** 3 # This is the correct function for volume mixing ratio
            if type(refractive_index).__name__ == 'DataFrame':
                refractive_index = refractive_index.iloc[:,0]
            # import pdb
            # pdb.set_trace()
            return out, n_mix(refractive_index, out)

    return out


def kappa_from_fofrh_and_sizedist(f_of_RH, dist, wavelength, RH, verbose = False, f_of_RH_collumn = None):
    """
    Calculates kappa from f of RH and a size distribution.
    Parameters
    ----------
    f_of_RH: TimeSeries
    dist: SizeDist
    wavelength: float
    RH: float
        Relative humidity at which the f_of_RH is taken.
    column: string
        when f_of_RH has more than one collumn name the one to be used
    verbose: bool

    Returns
    -------
    TimeSeries
    """

    def minimize_this(gf, sr, f_rh_soll, ext, wavelength, verbose = False):
        gf = float(gf)
        sr_g = sr.apply_growth(gf, how='shift_bins')
        sr_g_opt = sr_g.calculate_optical_properties(wavelength)
        ext_g = sr_g_opt.extinction_coeff_sum_along_d.data.values[0][0]
        f_RH = ext_g / ext
        out = f_RH - f_rh_soll

        if verbose:
            print('test growhtfactor: %s'%gf)
            print('f of RH soll/is: %s/%s'%(f_rh_soll,f_RH))
            print('diviation from soll: %s'%out)
            print('---------')
        return out

    # make sure f_of_RH has only one collumn
    if f_of_RH.data.shape[1] > 1:
        if not f_of_RH_collumn:
            txt = 'f_of_RH has multiple collumns (%s). Please name the one you want to use by setting the f_of_RH_collumn argument.'%(f_of_RH.data.columns)
            raise ValueError(txt)
        else:
            f_of_RH = f_of_RH._del_all_columns_but(f_of_RH_collumn)

    n_values = dist.data.shape[0]
    gf_calc = _np.zeros(n_values)
    kappa_calc = _np.zeros(n_values)
    f_of_RH_aligned = f_of_RH.align_to(dist)
    for e in range(n_values):
        frhsoll = f_of_RH_aligned.data.values[e][0]
        if _np.isnan(frhsoll):
            kappa_calc[e] = _np.nan
            gf_calc[e]  = _np.nan
            continue

        if type(dist.index_of_refraction).__name__ == 'float':
            ior = dist.index_of_refraction
        else:
            ior = dist.index_of_refraction.iloc[e][0]
        if _np.isnan(ior):
            kappa_calc[e] = _np.nan
            gf_calc[e]  = _np.nan
            continue

        sr = dist.copy()
        sr.data = sr.data.iloc[[e],:]
        sr.index_of_refraction = ior
        sr_opt = sr.calculate_optical_properties(wavelength)
        ext = sr_opt.extinction_coeff_sum_along_d.data.values[0][0]


        if ext == 0:
            kappa_calc[e] = _np.nan
            gf_calc[e]  = _np.nan
            continue

        if verbose:
            print('goal for f_rh: %s'%frhsoll)
            print('=======')

        gf_out = _fsolve(minimize_this, 1, args = (sr, frhsoll, ext, wavelength, verbose), factor=0.5, xtol = 0.005)
        gf_calc[e] = gf_out

        if verbose:
            print('resulting gf: %s'%gf_out)
            print('=======\n')

        kappa_calc[e] = kappa_simple(gf_out, RH, inverse=True)
    ts_kappa = _timeseries.TimeSeries(_pd.DataFrame(kappa_calc, index = f_of_RH_aligned.data.index, columns= ['kappa']))
    ts_kappa._data_period = f_of_RH_aligned._data_period
    ts_kappa._y_label = '$\kappa$'

    ts_gf = _timeseries.TimeSeries(_pd.DataFrame(gf_calc, index = f_of_RH_aligned.data.index, columns= ['growth factor']))
    ts_kappa._data_period = f_of_RH_aligned._data_period
    ts_gf._y_label = 'growth factor$'
    return ts_kappa, ts_gf