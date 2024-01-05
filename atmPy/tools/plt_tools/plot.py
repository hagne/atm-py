import numpy as _np
import matplotlib.pylab as _plt
import matplotlib.collections as _mcoll
import matplotlib.path as _mpath
from matplotlib import dates as _mdates
import pandas as _pd
from matplotlib.backend_bases import GraphicsContextBase as _GraphicsContextBase
from matplotlib.backend_bases import RendererBase as _RendererBase
import types as _types

#########
# The following hack makes the solid_joinstyle to round for linecollections ... https://stackoverflow.com/questions/11578760/matplotlib-control-capstyle-of-line-collection-large-number-of-lines
class GC(_GraphicsContextBase):
    def __init__(self):
        super().__init__()
        self._capstyle = 'round'

def custom_new_gc(self):
    return GC()

_RendererBase.new_gc = _types.MethodType(custom_new_gc, _RendererBase)
##########

def plot_gradiant_color(x, y, z = None, resample = 1, ax = None, colorbar = True, rescale= True, lc_kwargs = None):
    """ Plots a line given by x and y where the color changes along the line either gradually or according
    to z.

    Note
    ----
    based on following post:
    https://stackoverflow.com/questions/8500700/how-to-plot-a-gradient-color-line-in-matplotlib?answertab=active#tab-top

    Parameters
    ----------
    x: ndarray
    y: ndarray
    z:
    resample: int
        When there are too few data points the color resolution might be ugly. Resample will change the resolution.
    ax: matplotlib.axes._subplots.AxesSubplot instance
        If plotting on an existing axes is desired.
    rescale: bool
        if to recalculate limits of axes
    colorbar: bool
        If a colorbar is supposed to be plottet
    lc_kwargs: dict
        Dictionary containing kwargs of matplotlib.collections.LineCollection, e.g. linewidth, cmap, etc.

    Returns
    -------
    matplotlib.axes._subplots.AxesSubplot instance
    matplotlib.collections.LineCollection instance
    matplotlib.colorbar.Colorbar instance or None

    """


    def make_segments(x, y):
        """
        Create list of line segments from x and y coordinates, in the correct format
        for LineCollection: an array of the form numlines x (points per line) x 2 (x
        and y) array
        """

        points = _np.array([x, y]).T.reshape(-1, 1, 2)
        segments = _np.concatenate([points[:-1], points[1:]], axis=1)
        return segments

    if type(x[0]).__name__ == 'datetime64':
        x = _mdates.date2num(_pd.DatetimeIndex.to_pydatetime(_pd.to_datetime(x)))

    if type(y[0]).__name__ == 'datetime64':
        y = _mdates.date2num(_pd.DatetimeIndex.to_pydatetime(_pd.to_datetime(y)))

    if not ax:
        fig, ax = _plt.subplots()
    else:
        fig = ax.get_figure()

    if not lc_kwargs:
        lc_kwargs = {}

    if resample != 1:
        path = _mpath.Path(_np.column_stack([x, y]))
        verts = path.interpolated(steps=resample).vertices
        xv, yv = verts[:, 0], verts[:, 1]
    else:
        xv, yv = x, y

    # Default colors equally spaced on [0,1]:
    if z is None:
        z = _np.linspace(0.0, 1.0, len(xv))

    if type(z[0]).__name__ == 'datetime64':
        z = _mdates.date2num(_pd.DatetimeIndex.to_pydatetime(_pd.to_datetime(z)))

    # if resample is not 1 z needs to be resampled
    else:
        if resample != 1:
            sp = _np.linspace(0, 1, len(x))
            sp_new = _np.linspace(0, 1, len(xv))
            z = _np.interp(sp_new, sp, z)

    # Special case if a single number:
    # if not hasattr(z, "__iter__"):  # to check for numerical input -- this is a hack
    #     z = _np.array([z])

    # z = _np.asarray(z)

    segments = make_segments(xv, yv)
    lc = _mcoll.LineCollection(segments,
                               # colors = z, # for this to work colors would need to be generated first ... clim would not work anymore
                               array=z,
                               # norm=_plt.Normalize(0.0, 1.0),
                               **lc_kwargs)

    ax.add_collection(lc)
    xnn = xv[~ _np.isnan(xv)]
    ynn = yv[~ _np.isnan(yv)]
    if rescale:
        ax.set_ylim((ynn.min(), ynn.max()))
        ax.set_xlim((xnn.min(), xnn.max()))

    if colorbar:
        colorbar = fig.colorbar(lc)

    return ax ,lc ,colorbar