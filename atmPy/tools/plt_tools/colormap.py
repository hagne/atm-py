import numpy as _np
import matplotlib.colors as _mcolors
# import LinearSegmentedColormap as _LinearSegmentedColormap

import plt_tools._my_colormaps as my_colormaps

def creat_cmap(colors = None, color_range = 1, norm='linear', log_min=0.1, reverse=False):
    """
    Creates a custom colormap

    Parameters
    ----------
    colors: list of colors, default None will pick my favorite
    color_range: int
        when copying the color from somewhere else it is usually between 0 and 1 or 255. Set the upper end here.
    norm: string (['linear'], 'log')
        If log log_min has to be set.
    log_min: float (0,1)
        Below log_min there is one more color step all the way to 0
    reverse: bool
        If the colors aught to be reversed

    Returns
    -------
        matplotlib.colors.LinearSegmentedColormap instance

    """
    if type(colors).__name__ == 'NoneType':
        colors = _np.array([_np.array([0.0, 4., 76.]) / 255.,
                            _np.array([49., 130., 0.0]) / 255.,
                            _np.array([255, 197., 98.]) / 255.,
                            _np.array([245., 179., 223.]) / 255.,
                            _np.array([1, 1, 1]),
                            ])

    colors /= color_range

    if norm == 'linear':
        steps = _np.linspace(0, 1, len(colors))
    elif norm == 'log':
        if not 0 < log_min < 1:
            raise ValueError(' 0 < log_min < 1 has to be True')
        steps = _np.logspace(_np.log10(log_min), 0, len(colors))
        steps[0] = 0

    else:
        raise ValueError('norm has to be one of the following values: linear, log')

    if reverse:
        colors = colors[::-1]

    r = _np.zeros((len(colors), 3))
    r[:, 0] = steps
    r[:, 1] = colors[:, 0]
    r[:, 2] = colors[:, 0]

    g = _np.zeros((len(colors), 3))
    g[:, 0] = steps
    g[:, 1] = colors[:, 1]
    g[:, 2] = colors[:, 1]

    b = _np.zeros((len(colors), 3))
    b[:, 0] = steps
    b[:, 1] = colors[:, 2]
    b[:, 2] = colors[:, 2]

    cdict = {'red': r,
             'green': g,
             'blue': b
             }

    hag_cmap = _mcolors.LinearSegmentedColormap('hag_cmap', cdict)
    hag_cmap.set_bad('black')
    return hag_cmap


class FermiNormalize(_mcolors.Normalize):
    """creates a normalization for a colormap that is linear but with a step in it. It is the combination of an inverted fermi function and a linear function
    Usage
    -----
    ferminorm = FermiNormalize(vmin=0, vcenter=0.92, vmax=8, stepsteepness = 0.15, global_slope = 0.1)
    pc = a.pcolor(X,Y,Z, norm = ferminorm)
    """
    def __init__(self, vmin=None, vmax=None, vcenter=None, stepsteepness = 0.3, global_slope = 0.3, clip=False):
        self.vcenter = vcenter
        self.stepsteepness = stepsteepness
        _mcolors.Normalize.__init__(self, vmin, vmax, clip)

    def __call__(self, value, clip=None):
        # I'm ignoring masked values and all kinds of edge cases to make a
        # simple example...
#         x, y = [self.vmin, self.vcenter, self.vmax], [0, 0.25, 1]
        funk = lambda x : (1 / (_np.exp((self.vcenter-x)/self.stepsteepness) + 1)) + (0.03 * x)
        trans = (_np.ma.masked_array(funk(value)) - funk(self.vmin)) / (funk(self.vmax) - funk(self.vmin))
        return trans