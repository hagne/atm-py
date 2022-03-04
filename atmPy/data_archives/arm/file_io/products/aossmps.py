import numpy as _np
import xarray as _xr
from atmPy.data_archives.arm import _tools
from atmPy.aerosols.size_distribution import sizedistribution as _sd
from atmPy.general import timeseries as _ts
import atmPy
def open_path(path, window = None, average = None, verbose = True):
    """

    Parameters
    ----------
    path
    start_time
    end_time
    average: tuple [None]
        The purpose of this is to keep the memory usage low in case a lower reolution is required. E.g. (60, 's')

    Returns
    -------

    """
    def read_aossmps(file, verbose = False):
        ds = _xr.open_dataset(file, autoclose=True)
        data_dist = ds.number_size_distribution.to_pandas()
#         data_dist = data_dist.iloc[: ,:-1]
        # bincenters =
        bincenters = data_dist.columns.values
        #     dist = sd.SizeDist_TS(data_dist, bincenters, 'numberConcentration')
#         binedges = _np.unique(ds.aerodynamic_diameter_bound.data.flatten())[1:] * 1000
        binedges, bn = atmPy.aerosols.size_distribution.diameter_binning.bincenters2binsANDnames(bincenters)
        # normalize to sample flow rate
#         sample_flow_rate_cc_s = (ds.total_flow_rate.to_pandas() - ds.sheath_flow_rate.to_pandas()) * 1000 /60
#         data_dist = data_dist.divide(sample_flow_rate_cc_s, axis = 'index')
        out = {}
        out['data_dist'] = data_dist
        out['bincenters'] = bincenters
        out['binedges'] = binedges
        if verbose:
            print(file)
            print('shapes: {}, {}'.format(data_dist.shape, bincenters.shape))
        return out

    # start_time, end_time  = window
    files = _tools.path2filelist(path=path, window = window, product='aossmps')
    if verbose:
        print('Opening {} files.'.format(len(files)))
        print(_tools.path2info(files[0]))
    data_dist = None
    binedges = None
    for file in files:
        out = read_aossmps(file)
        ddt = _ts.TimeSeries(out['data_dist'])
        if average:
            ddt = ddt.average_time(average)
        ddt = ddt.data
        if isinstance(data_dist, type(None)):
            data_dist = ddt
            binedges = out['binedges']
        else:
            data_dist = data_dist.append(ddt, sort=True)
            # make sure bincenters did not change
            assert (_np.all(_np.equal(binedges, out['binedges'])))

    dist = _sd.SizeDist_TS(data_dist, binedges, 'dNdlogDp', ignore_data_gap_error=True)
    return dist