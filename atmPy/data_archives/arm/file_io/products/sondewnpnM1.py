import xarray as _xr
import atmPy.atmosphere.sounding as _sounding
import pandas as _pd

def read_netCDF(fname):
    # fname = '/Volumes/HTelg_4TB_Backup/arm_data/OLI/balloon_soundings/olisondewnpnM1.b1.20170401.172800.cdf'

    if type(fname) == str:
        fname = [fname]

    dfs = []
    for fn in fname:
        data = _xr.open_dataset(fn)
        df = _pd.DataFrame(data.alt.to_pandas(), columns=['Altitude'])

        df['Lat'] = data.lat.to_pandas()
        df['Lon'] = data.lon.to_pandas()

        df['Press'] = data.pres.to_pandas()
        df['Lon'] = data.lon.to_pandas()

        df['Dp'] = data.dp.to_pandas()
        df['Temp'] = data.tdry.to_pandas()
        df['Wind_speed'] = data.wspd.to_pandas()
        df['Wind_direction'] = data.deg.to_pandas()
        df['Wind_u'] = data.u_wind.to_pandas()
        df['Wind_v'] = data.v_wind.to_pandas()
        df['Wind_status'] = data.wstat.to_pandas()
        df['RH'] = data.rh.to_pandas()
        df['Ascent_rate'] = data.asc.to_pandas()

        dfs.append(df)

    df = _pd.concat(dfs).sort_index()
    out = _sounding.BalloonSounding(df)
    return out