from atmPy.general import measurement_site as _measurement_site
import pandas as _pd
import numpy as _np
from atmPy.radiation import solar as _solar
from atmPy.general import timeseries as _timeseries


class AOD_AOT(object):
    def __init__(self,
                 wavelengths = None,
                 site = None,
                 lat = None,
                 lon = None,
                 elevation = 0,
                 name = None,
                 name_short = None,
                 timezone = 0,
                 site_info = None):
        """This class is for column AOD or AOT meausrements at a fixed site. This class is most usfull for aerosol
        optical properties from a CIMEL (AERONET) or a MFRSR (SURFRAD)

        Parameters
        ----------
        wavelengths: dict
            Column names are often not reflecting the precise wavelength in the channel, but the typical wavelength.
            The dictionary translates column names to exact wavelength. If AOD is calculated and wavelengths is set
            wavelengths from this will be used instead of column names.
        site: atmPy.general.station instance
        lat, lon: location of site
            lat: deg north, lon: deg east
            e.g. lat = 40.05192, lon = -88.37309 for Bondville, IL, USA
        elevation: float
            elevation of site in meter (not very importent parameter, keeping it at zero is usually fine ...
            even vor Boulder
        name, name_short: str
            name and abbriviation of name for site
        timezone: int
            Timezon in houres, e.g. Bondville has the timezone -6.
            Note, this can lead to some confusion if the data is in UTC not in local time.... this needs to be improved
            ad some point"""

        self._aot = None
        self._aod = None
        self._sunposition = None
        self._timezone = timezone
        self.wavelengths = wavelengths

        if not isinstance(site, type(None)):
            self.site = site
        elif not isinstance(lat, type(None)):
            self.site = _measurement_site.Station(lat, lon, elevation, name=name, abbreviation=name_short, info = site_info)




    @property
    def sun_position(self):
        if not self._sunposition:
            if self._timezone != 0:
                date = self._timestamp_index +  _pd.to_timedelta(-1 * self._timezone, 'h')
            else:
                date = self._timestamp_index
            self._sunposition = _solar.get_sun_position(self.site.lat, self.site.lon, date)
            self._sunposition.index = self._timestamp_index
            self._sunposition = _timeseries.TimeSeries(self._sunposition)
        return self._sunposition

    @property
    def AOT(self):
        if not self._aot:
            if not self._aod:
                raise AttributeError('Make sure either AOD or AOT is set.')
            aot = self.AOD.data.mul(self.sun_position.data.airmass, axis='rows')
            aot.columns.name = 'AOT@wavelength(nm)'
            aot = _timeseries.TimeSeries(aot)
            self._aot = aot
        return self._aot

    @ AOT.setter
    def AOT(self,value):
        self._aot = value
        self._aot.data.columns.name = 'AOT@wavelength(nm)'
        self._timestamp_index = self._aot.data.index

    @property
    def AOD(self):
        if not self._aod:
            if not self._aot:
                raise AttributeError('Make sure either AOD or AOT is set.')
            aod = self.AOT.data.div(self.sun_position.data.airmass, axis='rows')
            aod.columns.name = 'AOD@wavelength(nm)'
            aod = _timeseries.TimeSeries(aod)
            self._aod = aod
        return self._aod

    @ AOD.setter
    def AOD(self,value):
        self._aod = value
        self._aod.data.columns.name = 'AOD@wavelength(nm)'
        self._timestamp_index = self._aod.data.index

    @property
    def ang_exp(self):
        return self._ang_exp

    @ang_exp.setter
    def ang_exp(self, value):
        self._ang_exp = value

    def aod2angstrom_exponent(self, column_1=500, column_2=870,
                              use_wavelength_from_column_names = None,
                              # wavelength_1=None, wavelength_2=None
                              ):
        """
        Calculates the angstrom exponents based on the AOD data.

        Parameters
        ----------
        column_1: type of column name
            column name of one of the two points used for the AOD calculation
        column_2: type of column name
            column name of the other of the two points used for the AOD calculation
        use_wavelength_from_column_names: bool [None]
            When the wavelength dictionary is set. Wavelengths from the dictionary are used instead of column names.
            Set this kwarg to True to ignore the wavelengths dictionary and use column names instead.

        Parameters (deprecated)
        -----------------------
        wavelength_1: float
            if the column name of column_1 is not accurate enough set the wavelenth used to calculate AOD here.
        wavelength_2: float
            as above for column_2

        Returns
        -------

        """
        if isinstance(self.wavelengths, type(None)) or use_wavelength_from_column_names:
            # if wavelength_1 == None:
            wavelength_1 = column_1
            # if wavelength_2 == None:
            wavelength_2 = column_2
        else:
            wavelength_1 = self.wavelengths[column_1]
            wavelength_2 = self.wavelengths[column_2]
        c1 = column_1
        c2 = column_2
        c1ex = wavelength_1
        c2ex = wavelength_2
        out = - _np.log10(self.AOD.data.loc[:, c1] / self.AOD.data.loc[:, c2]) / _np.log10(c1ex / c2ex)
        out = _timeseries.TimeSeries(_pd.DataFrame(out))
        setattr(self, 'ang_exp_{}_{}'.format(column_1, column_2), out)
        return out

    def cal_deviations(arglist):
        #     print(arglist)
        [idx, aod], [idxt, wls] = arglist

        #         wls = sfr.wavelength_matching_table.data.loc[idx,:]
        assert (len(wls) == len(aod))

        sfr.ang_exp.data.loc[idx, :]

        fitres = sp.stats.linregress(np.log10(wls.loc[[500, 870]]), np.log10(aod.loc[[500, 870]]))

        ff = lambda x: 10 ** (np.log10(x) * fitres.slope + fitres.intercept)

        aod_fit = ff(wls)
        diff = (aod - aod_fit)
        return diff

    def calculate_deviation_from_angstrom(sfr, wl1=500, wl2=870):
        """
        calcu
        Parameters
        ----------
        wl1
        wl2

        Returns
        -------

        """
        pool = mp.Pool(6)

        # out = pool.map(cal_deviations, sfr.AOD.data.iloc[135:140].iterrows())
        out = pool.map(cal_deviations, zip(sfr.AOD.data.iterrows(), sfr.wavelength_matching_table.data.iterrows()))

        # make dataframe from output
        df = pd.concat(out, axis=1)
        df = df.transpose()

        # xarray does not like timezone ... remove it
        df.index = pd.to_datetime(df.index.astype(str).str[:-6])

        # xarray does not like integer column ames
        df.columns = df.columns.astype(str)

        # other stuff
        df.sort_index(inplace=True)
        df.index.name = 'datetime'
        return df

