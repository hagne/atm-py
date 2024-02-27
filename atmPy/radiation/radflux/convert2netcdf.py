import pandas as pd
import xarray as xr
import numpy as np
import pathlib as pl

# Variable descriptions.
# Descriptions are mostly based on Chucks original readme document at
# https://gml.noaa.gov/aftp/data/radiation/surfrad/RadFlux/RadFlux_ReadMe.txt
variable_info = [
    {
        "old": "AU",
        "new": "distance_from_sun",
        "long_name": "earth-sun distance in AUs"
    },
    {
        "old": "SWdn",
        "new": "shortwave_down_best_estimate",
        "long_name": "best estimate downwelling SW from sum or global pyranometer (W/m^2)"
    },
    {
        "old": "CSWdn",
        "new": "shortwave_down_clear_sky_estimate",
        "long_name": "estimated clear-sky downwelling SW (W/m^2)"
    },
    {
        "old": "LWdn",
        "new": "longwave_down",
        "long_name": "downwelling LW from pyrgeometer (W/m^2)"
    },
    {
        "old": "CLWdn",
        "new": "longwave_down_clear_sky_estimate",
        "long_name": "estimated clear-sky downwelling LW (W/m^2)"
    },
    {
        "old": "SWup",
        "new": "shortwave_up",
        "long_name": "upwelling SW from pyranometer (W/m^2)"
    },
    {
        "old": "CSWup",
        "new": "shortwave_up_clear_sky_estimate",
        "long_name": "estimated clear-sky upwelling SW (W/m^2)"
    },
    {
        "old": "LWup",
        "new": "longwave_up",
        "long_name": "upwelling LW from pyrgeometer (W/m^2)"
    },
    {
        "old": "CLWup",
        "new": "longwave_up_clear_sky_estimate",
        "long_name": "estimated clear-sky upwelling LW (W/m^2)"
    },
    {
        "old": "DifSW",
        "new": "shortwave_diffuse_down_best_estimate",
        "long_name": "measured downwelling diffuse SW (W/m^2)"
    },
    {
        "old": "CDifSW",
        "new": "shortwave_diffuse_down_clear_sky_estimate",
        "long_name": "estimated clear-sky downwelling diffuse SW (W/m^2)"
    },
    {
        "old": "DirSW",
        "new": "shortwave_direct",
        "long_name": "measured downwelling direct SW (W/m^2)"
    },
    {
        "old": "CDirSW",
        "new": "shortwave_direct_clear_sky_estimate",
        "long_name": "estimated clear-sky downwelling direct SW (W/m^2)"
    },
    {
        "old": "ClrF",
        "new": "flag_clear_sky",
        "long_name": "Clear sky flag, 1 if SW detected clear sky, 2 if LW detected, 9 if CLW>LW, 3 if only std and Ta-Te diff OK and ONLY LWup accepted as clear LWup [NOT LWdn!!!], else 0 if cloudy"
    },
    {
        "old": "TauF",
        "new": "flag_tau",
        "long_name": "Tau flag, 1 if liq g used, 2 if ice g used, 0 if not calculated"
    },
    {
        "old": "TlmF",
        "new": "flag_tlim",
        "long_name": "T limit flag, 1 if SW Scv used, 2 if LW Scv used, 3 if avg Ec used, 4 if lim=0.965*Ta used, 5 if just config limit temp used, 0 if not calculated"
    },
    {
        "old": "LWScv",
        "new": "sky_cover_longwave",
        "long_name": "estimated effective LW fractional sky cover"
    },
    {
        "old": "SWScv",
        "new": "sky_cover_shortwave",
        "long_name": "estimated fractional sky cover from SW"
    },
    {
        "old": "CldTau",
        "new": "cloud_OD",
        "long_name": "estimated effective visible cloud optical depth  (only for SWScv>0.95)"
    },
    {
        "old": "CldTrn",
        "new": "cloud_transmissivity",
        "long_name": "estimated effective SW cloud transmissivity (SWdn/CSWdn ratio)"
    },
    {
        "old": "TeLim",
        "new": "cloud_ice_temp_limit",
        "long_name": "Ice cloud temp limit (K)"
    },
    {
        "old": "LWTe",
        "new": "sky_brightness_temp",
        "long_name": "Sky brightness temp from LWdn (K)"
    },
    {
        "old": "CldTmp",
        "new": "cloud_radiating_temp",
        "long_name": "estimated effective cloud radiating temperature"
    },
    {
        "old": "CldHgt",
        "new": "cloud_radiating_height",
        "long_name": "estimated effective cloud radiating height"
    },
    {
        "old": "Tair",
        "new": "temp_air",
        "long_name": "air temperature (K)"
    },
    {
        "old": "VPrs",
        "new": "vapor_pressure",
        "long_name": "vapor pressure  (mb)"
    },
    {
        "old": "RH",
        "new": "rh",
        "long_name": "Relative Humidity (%)"
    },
    {
        "old": "RHfac",
        "new": "rh_frac",
        "long_name": "RH adjustment to Ec"
    },
    {
        "old": "Ec",
        "new": "longwave_emissivity_clear_sky",
        "long_name": "effective clear-sky LW emissivity"
    },
    {
        "old": "Wspd",
        "new": "wind_speed",
        "long_name": "Wind speed (same as input)"
    },
    {
        "old": "LWlw",
        "new": "correction_term_lwup_lwdn",
        "long_name": "(if included) Contribution to clear-sky LWup from LWdn term (W/m^2)"
    },
    {
        "old": "SWlw",
        "new": "correction_term_lwup_swnet",
        "long_name": "(if included) Contribution to clear-sky LWup from SWnet term (W/m^2)"
    },
    {
        "old": "RHlw",
        "new": "correction_term_lwup_rh",
        "long_name": "(if included) Contribution to clear-sky LWup from RH term (W/m^2)"
    },
    {
        "old": "Wslw",
        "new": "correction_term_lwup_wspd",
        "long_name": "(if included) Contribution to clear-sky LWup from Wspd term (W/m^2)"
    },
    {
        "old": "aprs",
        "new": "aprs",
        "long_name": "N/A"
    },
    {
        "old": "uvbdn",
        "new": "uvbdn",
        "long_name": "N/A"
    },
    {
        "old": "pardn",
        "new": "pardn",
        "long_name": "N/A"
    },
    {
        "old": "wdir",
        "new": "wdir",
        "long_name": "N/A"
    },
    {
        "old": "gcorr",
        "new": "gcorr",
        "long_name": "N/A"
    },
    {
        "old": "lwdtc",
        "new": "lwdtc",
        "long_name": "N/A"
    },
    {
        "old": "lwdtd",
        "new": "lwdtd",
        "long_name": "N/A"
    },
    {
        "old": "lwutc",
        "new": "lwutc",
        "long_name": "N/A"
    },
    {
        "old": "lwutd",
        "new": "lwutd",
        "long_name": "N/A"
    },
    {
        "old": "CosZ",
        "new": "solar_zenith_angle",
        "long_name": "Angle between Sun and the vertical",
        "unit": "degree"
    },
    {
        "old": "time_local",
        "new": "time_local",
        "long_name": "Local standard time"
    }
]


def chuck_file2dataset(path2file, verbose=False):
    p2f = path2file
    
    df = pd.read_csv(p2f, delim_whitespace=True)

    # at some point the column names went from upper and lower case to just lower case
    df.columns = [c.lower() for c in df.columns]
    
    df.index = df.apply(lambda row: pd.to_datetime(f'{row.zdate:0.0f}-{row.ztim:04.0f}'), 
                        axis=1)
    df.index.name = 'time'
    
    df['time_local'] = df.apply(lambda row: pd.to_datetime(f'{row.ldate:0.0f}-{row.ltim:04.0f}'), axis=1)
    
    df.drop(['zdate', 'ztim', 'ldate', 'ltim'], axis=1, inplace=True)
    
    # instead of giving the cosine of the zenith angle I rather give the zenith angle
    df.cosz = np.rad2deg(np.arccos(df.cosz))
    
    # lets set invalid data to nan
    df[df == -9999] = np.nan
    
    ds = xr.Dataset()
    cols = list(df.columns)
    optional = ['lwlw', 'swlw', 'rhlw', 'wslw']
    for var in variable_info:
        # break
        try:
            ds[var['new']] = df[var['old'].lower()]
        except KeyError as e:
            if e.args[0] in optional:
                if verbose:
                    print(f'Optional variable {e.args[0]} not found')
                # cols.pop(cols.index(var['old'].lower()))
                continue
            else:
                raise
        for a in var:
            if a in ['new', 'old']:
                continue
            ds[var['new']].attrs[a] = var[a]
        cols.pop(cols.index(var['old'].lower()))

    assert (len(cols) == 0), f'Not all columns where assigned. Collumns left over: {cols}.'

    return ds


class Convert(object):
    def __init__(self,
                 path2fld_in='/nfs/iftp/aftp/g-rad/surfrad/RadFlux/',
                 path2fld_out='/nfs/grad/surfrad/products_level4/radflux/v{version}/',
                 sites=['tbl', 'dra', 'fpk', 'gwn', 'psu', 'sxf', 'bon']
                 ):
        self.version = '1.0'
        self.p2fld_in = pl.Path(path2fld_in)
        self.p2fld_out = pl.Path(path2fld_out.format(version=self.version))
        self.sites = sites
        self._workplan = None

    def process(self, verbose=False, error_handling='raise'):
        self.p2fld_out.mkdir(exist_ok=True, parents=True)
        for idx, row in self.workplan.iterrows():
            try:
                if verbose:
                    print(row.p2f_in)
                ds = chuck_file2dataset(row.p2f_in)
                row.p2f_out.parent.mkdir(exist_ok=True, parents=True)
                ds.to_netcdf(row.p2f_out)
            except: #Exception as e:
                if error_handling == 'skip':
                    continue
                else:
                    print(f'proplem with {row.p2f_in}.')
                    raise

    @property
    def workplan(self):
        if isinstance(self._workplan, type(None)):
            wp = pd.DataFrame()
            for fld in self.p2fld_in.glob('*'):
                if not fld.is_dir():
                    continue
                elif fld.name not in self.sites:
                    continue
                wpt = pd.DataFrame(fld.glob('**/*.lw1'), columns=['p2f_in'])
                wpt['site'] = fld.name
                wp = pd.concat([wp, wpt])
            
            # generate the time (or date) index
            wp.index = wp.apply(lambda row: pd.to_datetime(row.p2f_in.name[:8]), 
                                axis=1)
            wp.index.name = 'date'
            
            # row = wp.iloc[0]
            
            # generate the output paths
            wp['p2f_out'] = wp.apply(
                lambda row:
                    self.p2fld_out
                    .joinpath(row.site)
                    .joinpath(str(row.name.year))
                    .joinpath(f'{row.site}.radflux.v{self.version}.{row.p2f_in.name[:8]}.nc'),
                axis=1
            )
            
            # remove if output file exists
            wp = wp[~(wp.apply(lambda row: row.p2f_out.is_file(), axis=1))]
            wp.sort_index(inplace = True)
            self._workplan = wp
        return self._workplan

    @workplan.setter
    def workplan(self, value):
        self._workplan = value





















