from mpl_toolkits.axes_grid1.inset_locator import inset_axes as _inset_axes
from mpl_toolkits.axes_grid1 import make_axes_locatable as _make_axes_locatable



def colorbar_axis_split_off(mappable, ax, position = 'right', size = '5%', pad = 0.1, cb_kwargs = None):
    """
    Splits of a section of the axis and uses it as the colorbar. This way the colorbar will have the same hight or width
    as the original axis. It also makes it easier to sey where the colorbar goes if the figure contains multiple axes.

    Parameters
    ----------
    mappable: mappable instance or None
        e.g. what is returned by imshow() or pcolormesh(). If None a pice of axes will still be cut off but no colorbar
        is plotted. This is usefull if you share axes and you want to apply to change to all axes.
    ax: AxesSubplot
    position: string (['right'], 'left', 'bottom', 'top')
        Where to split off the section of the axis
    size: axes_grid.axes_size compatible
        e.g. '5%'
    pad: float
        distance of split off from rest of axis
    cb_kwargs: dict [None]
        pass colorbar key word arguments in a dictionary here

    Returns
    -------
    matplotlib.colorbar.Colorbar
    matplotlib.axes._axes.Axes
        This is the axis which was split off and which contains the colorbar

    """
    f = ax.get_figure()
    divider = _make_axes_locatable(ax)
    cax = divider.append_axes(position, size=size, pad=pad)
    if isinstance(mappable, type(None)):
        cax.axis('off')
        return None, cax
    if not cb_kwargs:
        cb_kwargs = {}
    if position in ['top', 'bottom']:
        if 'orientation' not in cb_kwargs.keys():
            cb_kwargs['orientation'] = 'horizontal'
    cb = f.colorbar(mappable, cax=cax, **cb_kwargs)
    return cb, cax

def colorbar_inside_plot(ax, mappable, extend = ('17%', '60%'), extend_cb = ('25%', '90%'), loc = 7, loc_cb = 6, color_bg = [1, 1, 1, 0.5], colorbar_kw = {}):
    """
    Plots a colorbar inside the current plot.
    Parameters
    ----------
    ax: AxesSubplot
    mappable: mappable instance
        e.g. what is returned by imshow() or pcolormesh()
    extend: tuple
        The extend the background of the colorbar is taking up of the underlying axis in %
    extend_cb: tuple
        The extend the colorbar is taking up inside the background axis
    loc: int
        Location option as for legend see https://stackoverflow.com/questions/28521744/error-of-adding-a-legend-for-a-plot-in-python-3-2-matplotlib
    loc_cb: int
        Location of colorbar inside background axis
    color_bg: what ever matplotlib excepts for colors
        Color of background.

    Returns
    -------
    Colorbar instance
    AxesHostAxes: Axis that containst the colorbar
    AxesHostAxes: Axis that makes up the background
    """

    # first insert that acts as a background of the colorbar ... so you can see the text on the underlying plot
    cbbox = _inset_axes(ax, extend[0], extend[1], loc = loc)
    [cbbox.spines[k].set_visible(False) for k in cbbox.spines]  # remove axes
    cbbox.tick_params(axis='both',
                      left='off',
                      top='off',
                      right='off',
                      bottom='off',
                      labelleft='off',
                      labeltop='off',
                      labelright='off',
                      labelbottom='off')                        # remove ticks
    
    cbbox.set_facecolor(color_bg)                            # make transparent
    cbbox.grid(False)
    cbbox.set_axis_off() # there where still some ticklabel showing up, even with the tickparams set above?!?
    # second insert that will be filled by the colorbar
    # extend_cb = ('25%', '90%')
    ax_cb = _inset_axes(cbbox, extend_cb[0], extend_cb[1], loc = loc_cb)
    f = ax.get_figure()
    cb=f.colorbar(mappable, cax=ax_cb, **colorbar_kw) #make colorbar
    return cb, ax_cb, cbbox

