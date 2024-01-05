

def set_shared_label(a, label, axis='x', labelpad=0.01):
    """Set a single shared label for multiple axes

    Parameters
    ----------
    a: list of matplotlib.axes._subplots.AxesSubplot
    label: string
    axis: string (['x'], 'y')
    labelpad: float
        Padding of label with respect to ticklabels

    Returns
    -------

    """


    f = a[0].get_figure()
    f.canvas.draw()  # sets f.canvas.renderer needed below

    if axis == 'y':
        # get the center position for all plots
        top = a[0].get_position().y1
        bottom = a[-1].get_position().y0

        # get the coordinates of the left side of the tick labels
        x0 = 1
        for at in a:
            at.set_ylabel('')  # just to make sure we don't and up with multiple labels
            bboxes, _ = at.yaxis.get_ticklabel_extents(f.canvas.renderer)
            # bboxes = bboxes.inverse_transformed(f.transFigure) #deprecated
            bboxes = bboxes.transformed(f.transFigure.inverted())
            xt = bboxes.x0
            if xt < x0:
                x0 = xt
        tick_label_left = x0

        # set position of label
        a[-1].set_ylabel(label)
        a[-1].yaxis.set_label_coords(tick_label_left - labelpad, (bottom + top) / 2, transform=f.transFigure)
        labelinstance = a[-1].get_ylabel()

    elif axis == 'x':
        right = a[0].get_position().x1
        left = a[-1].get_position().x0

        y0 = 1
        for at in a:
            at.set_xlabel('')  # just to make sure we don't and up with multiple labels
            bboxes, _ = at.xaxis.get_ticklabel_extents(f.canvas.renderer)
            bboxes = bboxes.inverse_transformed(f.transFigure)
            yt = bboxes.y0
            if yt < y0:
                y0 = yt
        tick_label_bottom = y0

        a[-1].set_xlabel(label)
        a[-1].xaxis.set_label_coords((right + left) / 2, tick_label_bottom - labelpad, transform=f.transFigure)  # tick_label_left - labelpad,(bottom + top)/2, transform=f.transFigure)
        labelinstance =  a[-1].get_xlabel()

    return labelinstance
