import xarray as _xr
import numpy as _np
import atmPy.aerosols.size_distribution.sizedistribution as _sizedist
import pandas as _pd

def read_netCDF(fname):
    # fname = '/Volumes/HTelg_4TB_Backup/arm_data/OLI/uhsas/oliaosuhsasM1.a1.20170401.000008.nc'

    if type(fname) == str:
        fname = [fname]

    sds = []
    for fn in fname:
        data = _xr.open_dataset(fn)

        if data.sampling_interval.split()[1] != 'seconds':
            raise ValueError('This should be seconds, but it is {}.'.format(data.sampling_interval.split()[1]))

        if not _np.all((data.upper_size_limit.data[:-1] - data.lower_size_limit.data[1:]) == 0):
            raise ValueError('Something is wrong with the bins')

        # flow rate variable name changed at some point
        if 'sample_flow_rate' in data.variables.keys():
            flowrate = data.sample_flow_rate
        elif 'sampling_volume' in data.variables.keys():
            flowrate = data.sampling_volume

        if flowrate.units not in ['sccm', 'cc/min']:
            raise ValueError('units has to be sccm, but is {}.'.format(flowrate.units))

        sd = data.size_distribution.to_pandas()

        # normalize total numbers to numbers/(cc)
        ## normalize to integration interval
        sd /= float(data.sampling_interval.split()[0])  # normalize to integration interval

        ## normalize to flow rate
        flowrate = flowrate.values / 60.
        sd = sd.divide(flowrate, axis=0)

        sds.append(sd)

    sd = _pd.concat(sds).sort_index()


    binedges = _np.append(data.lower_size_limit.data, data.upper_size_limit.data[-1])
    sdts = _sizedist.SizeDist_TS(sd, binedges, 'numberConcentration')
    sdts._data_period = float(data.sampling_interval.split()[0])
    sdts = sdts.convert2dNdlogDp()
    return sdts