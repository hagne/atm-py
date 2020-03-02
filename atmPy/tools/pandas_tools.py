import numpy as np
import matplotlib.pylab as plt
import pandas as _pd

def ensure_column_exists(df, col_name, col_alt = False, raise_error = True):
    """Checks if a particular name is among the column names of a dataframe. Alternative names can be given, which when
    found will be changed to the desired name. The DataFrame will be changed in place. If no matching name is found an
    AttributeError is raised.

    Parameters
    ----------
    col_name: pandas.DataFrame
    col_alt: bool or list
        list of aternative names to look for

    """
    out = True
    if col_name not in df.columns:
        renamed = False
        if col_alt:
            for dta in col_alt:
                if dta in df.columns:
                    df.rename(columns={dta:col_name}, inplace = True)
                    renamed = True
                else:
                    pass
        if not renamed:
            txt = 'Column %s not found.'%col_name
            if col_alt:
                txt += 'Neither one of the alternatives: %s'%(col_alt)
            out = False
            if raise_error:
                raise AttributeError(txt)
    return out

def plot_dataframe_meshgrid(df, xaxis = 0, ax = None, pc_kwargs = {}, cb_kwargs = {}):
    axes_list = [df.index, df.columns]
    x_index = axes_list[xaxis]

    if xaxis !=0:
        df = df.swapaxes(0,1)
        y_index = axes_list[0]
    else:
        y_index = axes_list[1]

    z = df.values.transpose()
    x = np.repeat(np.array([x_index]), y_index.shape[0], axis = 0)
    y = np.repeat(np.array([y_index]), x_index.shape[0], axis = 0).transpose()

    if ax:
        a = ax
        f = a.get_figure()
    else:
        f,a = plt.subplots()

    if not pc_kwargs:
        pc_kwargs = {}

    pc = a.pcolormesh(x, y , z, **pc_kwargs)


    if 'datetime' in df.index.dtype.__str__():
        f.autofmt_xdate()

    if cb_kwargs:
        # if passing an empty dict the colorbar is not plot ->
        if type(cb_kwargs) == bool or len(cb_kwargs) == 0:
            cb_kwargs = {'use_gridspec': True}
        cb = f.colorbar(pc, **cb_kwargs)
    else:
        cb = None

    a.set_xlabel(df.index.name)
    a.set_ylabel(df.columns.name)

    # nans, screw up the scaling, therefore ...
    # if np.any(np.isnan(df.values)):
    #     values = df.values
    #     values = values[~ np.isnan(values)]
    #     pc.set_clim((values.min(),values.max()))
    return f,a,pc,cb

def plot_panel_meshgrid(panel, xaxis = 0, yaxis = 1, sub_set = 0, ax = None, kwargs = {}):

    valid_axes = np.array([0,1,2])
    zaxis = valid_axes[np.logical_and(valid_axes != xaxis, valid_axes != yaxis)][0]
    axes_list = [panel.items, panel.major_axis, panel.minor_axis]

    axes_idx_list_tobe = np.array([int(xaxis),int(yaxis),int(zaxis)])
    axes_idx_list_is = valid_axes.copy()
    x_index = axes_list[xaxis]
    y_index = axes_list[yaxis]

    if axes_idx_list_tobe[0] != axes_idx_list_is[0]:
        axes_idx_list_is[0], axes_idx_list_is[axes_idx_list_tobe[0]] = axes_idx_list_is[axes_idx_list_tobe[0]], axes_idx_list_is[0]
        panel = panel.swapaxes(0,axes_idx_list_tobe[0])
    if axes_idx_list_tobe[1] != axes_idx_list_is[1]:
        axes_idx_list_is[1], axes_idx_list_is[2] = axes_idx_list_is[2], axes_idx_list_is[1]
        panel = panel.swapaxes(1,2)
    if not np.array_equal(axes_idx_list_is,axes_idx_list_tobe):
        txt = 'not possible'
        raise ValueError(txt)

    z = panel.values[:,:,sub_set].transpose()
    x = np.repeat(np.array([x_index]), y_index.shape[0], axis = 0)
    y = np.repeat(np.array([y_index]), x_index.shape[0], axis = 0).transpose()

    if ax:
        a = ax
        f = a.get_figure()
    else:
        f,a = plt.subplots()

    pc = a.pcolormesh(x, y , z, **kwargs)


    if 'datetime' in panel.items.dtype_str:
        f.autofmt_xdate()
    cb = f.colorbar(pc)
    a.set_xlabel(panel.items.name)
    a.set_ylabel(panel.major_axis.name)
    cb.set_label(panel.minor_axis[sub_set])
    pc.set_clim(z[~ np.isnan(z)].min(), z[~ np.isnan(z)].max())
    return f,a,pc,cb


def find_contour(df, value, xaxis = 0):
    """find a contour (line of equal values in a 2D array)"""
    axes_list = [df.index, df.columns]
    x_index = axes_list[xaxis]

    if xaxis !=0:
        df = df.swapaxes(0,1)
        y_index = axes_list[0]
    else:
        y_index = axes_list[1]

    z = df.values.transpose()
    x = np.repeat(np.array([x_index]), y_index.shape[0], axis = 0)
    y = np.repeat(np.array([y_index]), x_index.shape[0], axis = 0).transpose()


    timestamp = np.zeros(z.shape[1])
    timestamp.dtype = x.dtype
    idx = np.zeros(z.shape[1])
    dt = np.zeros(z.shape[1])
    for e, zt in enumerate(z.transpose()):
        temp = abs(zt[~np.isnan(zt)] - value)
        if temp.shape[0] == 0:
            dt[e] = np.nan
        else:
            idx[e] = temp.argmin()
            dt[e] = float(y[int(idx[e]),e])
        timestamp[e] = x[0,e]

    contour = _pd.DataFrame(dt, index = timestamp, columns = ['%s_contour'%value])
    return contour