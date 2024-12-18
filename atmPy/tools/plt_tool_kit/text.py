import numpy as _np
import pandas as _pd
import matplotlib.pyplot as plt

def latex_formula(formula, textkwargs = {}):
    """
    Generates a figure with only the latex formula, noting else!
    Example
    --------
    t = f'$\sigma = \sqrt{{\sigma_{{\\uparrow}}^2 + \sigma_{{\\downarrow}}^2}} = {sigupdown:0.1f}$'
    to = pltkit.text.latex_formula(t)
    
    """
    # Create a figure and an axes

    if 'fontsize' not in textkwargs:
        textkwargs={'fontsize':40}
    
    f, a = plt.subplots()

    # Add text at the center of the axes
    to = a.text(0, 0, formula, ha='center', va='center', **textkwargs)
    to.set_horizontalalignment('left')
    # Turn off the axes
    a.set_axis_off()

    # Draw the canvas to update the renderer and get the correct bounding box
    f.canvas.draw()

    # Get the bounding box of the text in inches
    bbox = to.get_window_extent().transformed(f.dpi_scale_trans.inverted())

    # Set the figure size based on the text bounding box, adding a small margin
    margin = 0.  # inches
    f.set_size_inches(bbox.width + 2 * margin, bbox.height + 2 * margin)

    # Adjust the layout to remove any extra whitespace
    f.tight_layout(pad=0)

    # Show or save the figure
    # plt.show()
    return to 



def get_angle_of_line_at_pos(data, x, aspectratio=1):
    """Returns the angle of the line at the given index"""
    px = x
    try:
        py = data.loc[x, 'y']
    except IndexError:
        raise IndexError('The postion is outside the date range ... adjust position to make this work')
    idx = data.index.get_loc(x)-1
    psx = data.iloc[idx].name
    psy = data.iloc[idx].y
    dx = px - psx
    dy = py - psy
    dy *= aspectratio
    deg = _np.degrees(_np.arctan2(dy , dx))
    return deg


def get_interpolated_datapoint(x, y, position, axes='x', return_data=False, verbose = True):
    data = _pd.DataFrame({'y': y}, index = x)
    if position not in data.index:
        data.loc[position] = _np.nan
        data.sort_index(inplace = True)
        data = data.interpolate(method = 'index')
    data.sort_index(inplace = True)
    px = position
    py = data.loc[position, 'y']
    if return_data:
        return px, py, data
    else:
        return px, py

def get_aspect(ax):
    fig = ax.figure
    ll, ur = ax.get_position() * fig.get_size_inches()
    width, height = ur - ll
    axes_ratio = height / width
    aspect = axes_ratio / ax.get_data_ratio()

    return aspect

def add_text_along_graph(graph, txt, position, axes='x', txtkwargs = None,  bbox = None):
    """as the name says :-)

    Parameters
    ----------
    graph: matplotlib.lines.Line2D
    txt: str
    position: float
        this is the position with respect to the axes devined in axes
    axes: 'x' or 'y'
        Which axes the position argument refers to.
    txtkwargs: dict
        kwargs passed when plt.text is called.
        Default: dict(va='center', ha='center')
        More info at:
        https://matplotlib.org/3.1.0/api/_as_gen/matplotlib.axes.Axes.text.html
    bbox: dict
        Arguments for a box around the text. If none desired set to False.
        Default: {'boxstyle': 'round,pad=0.2'}
        More info at:
        https://matplotlib.org/3.1.1/api/_as_gen/matplotlib. ches.BoxStyle.html

    Returns
    -------

    """
    a = graph._axes
    f = a.get_figure()
    col = graph.get_color()
    xt, yt = graph.get_data()

    bboxecef = False
    if isinstance(bbox, type(None)):
        bboxecef = True
        bbox = {'boxstyle': 'round,pad=0.2',
                'linewidth': 0.5,
                'facecolor': 'white',
                # 'edgecolor':
                }
    elif bbox == False:
        bbox = None
    else:
        assert(isinstance(bbox, dict))

    if isinstance(txtkwargs, type(None)):
        txtkwargs = dict(va='center', ha='center')
    else:
        assert(isinstance(txtkwargs, dict))


    px, py, data = get_interpolated_datapoint(xt, yt, position, axes, return_data=True)
    txo = a.text(px, py, txt, bbox=bbox, **txtkwargs)
    aspect_ratio = get_aspect(a)
    txo.set_rotation(get_angle_of_line_at_pos(data, px, aspectratio = aspect_ratio))
    boxo = txo.get_bbox_patch()
    if bboxecef:
        boxo.set_edgecolor(col)
    return txo, boxo