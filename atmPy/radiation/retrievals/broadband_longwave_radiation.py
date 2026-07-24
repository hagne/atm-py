
import xarray as xr
import numpy as np
import scipy
import atmPy.general.measurement_site as atmgms


default_config = dict()

class _DatasetRef(object):
    def __init__(self, dataset):
        """A simple class to hold a reference to a dataset. This is used to allow multiple classes to share the same dataset without having to pass it around explicitly."""
        self.dataset = dataset

class LongwaveIrradiation(object):
    def __init__(self, dataset, site = None, verbose = False, _dataset_ref = None):
        """Base class for longwave irradiation datasets. Provides shared functionality and properties for downwelling and upwelling longwave irradiation datasets.
        Parameters
        ----------
        dataset : xr.Dataset
            A dataset that contains the relevant irradiation variable(s) and the coordinates 'datetime'.
        site : atmPy.general.measurement_sites.Station, atmPy.general.measurement_sites.MovingPlatform, optional
            A measurement site object that can be used to calculate sun position and other site-specific parameters. 
            """
        self.verbose = verbose  
        self._dataset_ref = _dataset_ref if _dataset_ref is not None else _DatasetRef(dataset)
        if _dataset_ref is None:
            self.dataset = dataset
        self._site = site
        self._sun_position_variables = ['zenith', 'zenith_geometric', 'elevation_geometric', 'elevation', 'azimuth', 'equation_of_time', 'airmass', 'airmass_absolute', 'sun_earth_distance']
        self._sun_position_variables_ds = [f'solar_{v}' for v in self._sun_position_variables] # for internal use, to avoid name clashes with potential variables in the dataset.
        assert('datetime' in dataset or 'datetime' in dataset.dims), 'Time coordinate has to be called datetime .... sorry, i know that is an unconventional choise for the time cooridinate name.'

    @property
    def site(self):
        if isinstance(self._site, type(None)):
            """Todo: This is just copyied over from the shortwave irradiation class, centralize this in the future."""
            if self.verbose:
                print('No site provided, trying to infer site information from dataset.')
            if 'latitude' in self.dataset.attrs and 'longitude' in self.dataset.attrs:
                if self.verbose:
                    print('Found latitude and longitude in dataset attributes. Assuming fixed station.')
                if 'altitude' in self.dataset.attrs:
                    alt = self.dataset.attrs['altitude']
                else:
                    alt = 0
                self._site = atmgms.Station(lat = self.dataset.attrs['latitude'], lon = self.dataset.attrs['longitude'], alt = alt)
            elif 'latitude' in self.dataset and 'longitude' in self.dataset:
                if self.verbose:
                    print('Found latitude and longitude in dataset variables. ', end = '')
                if 'altitude' in self.dataset:
                    alt = self.dataset.altitude
                else:
                    alt = 0
                if 'datetime' in self.dataset.latitude.dims:
                    if self.verbose:
                        print('Latitude and longitude are time dependent. Assuming moving platform.')
                    self._site = atmgms.MovingPlatform(lat = self.dataset['latitude'], lon = self.dataset['longitude'], alt = alt)
                else:
                    if self.verbose:
                        print('Latitude and longitude are not time dependent. Assuming fixed station.')
                    self._site = atmgms.Station(lat = self.dataset['latitude'], lon = self.dataset['longitude'], alt = alt)

            else:
                raise ValueError('No site information found. Set site keyword or provide latitude and longitude in dataset attributes or variables. Note, time dependent latitude and longitude will be interpreted as a moving platform.')
        return self._site

            
            
    @property
    def dataset(self):
        return self._dataset_ref.dataset

    @dataset.setter
    def dataset(self, value):
        self._dataset_ref.dataset = value

    def drop_vars(self, names, *, errors = 'raise'):
        """Drop variables and update the shared dataset reference. In principle, this is an in place operation."""
        self.dataset = self.dataset.drop_vars(names, errors = errors)
        return self.dataset

    def get_attr(self, attr):
        if attr not in self.dataset.attrs:
            if self.verbose:
                print(f'{attr} attribute is not set, using default')
            self.dataset.attrs[attr] = default_config[attr]
        return self.dataset.attrs[attr]

    
    @property
    def sun_position(self):
        if not np.all([v in self.dataset for v in self._sun_position_variables_ds]):
        # if isinstance(self._sun_position, type(None)):
            sp = self.site.get_sun_position(self.dataset.datetime)
            self.tp_sp = sp
            for v in self._sun_position_variables:
                self.dataset[f'solar_{v}'] = sp[v]
        return self.dataset[self._sun_position_variables_ds]
    
    @sun_position.setter
    def sun_position(self, value):
        self.dataset.datetime.identical(value.datetime)
        for v in self._sun_position_variables:
            self.dataset[f'solar_{v}'] = value[v]

    @property
    def met_data(self):
        # if isinstance(self._metdata, type(None)):
        if 'pressure' not in self.dataset:
            raise ValueError('Metdata is not present in the dataset. Set it by passing data or filename.')
        return self.dataset[['pressure','temperature']]
    
    @met_data.setter
    def met_data(self, value):
        assert(False), 'Setting met data is not implemented yet.'


class DownwellingLongwaveIrradiation(LongwaveIrradiation):
    def __init__(self, dataset, *args, **kwargs):
        super().__init__(dataset, *args, **kwargs)
        if 'downwelling_longwave_irradiation' not in self.dataset:
            raise ValueError('downwelling_longwave_irradiation is not present in the dataset. Set it by passing data or filename.')
        



    @property
    def sky_brightness_temperature(self):
        if 'sky_brightness_temperature' not in self.dataset:
            self.dataset['sky_brightness_temperature'] = self._calculate_sky_brightness_temperature()
        return self.dataset['sky_brightness_temperature']
    

    def _calculate_sky_brightness_temperature(
        self,
    ) -> xr.Dataset:
        """Calculate effective sky brightness temperature and air–sky difference.
        References
        ----------
        Long, C. N., and D. D. Turner (2008).
        A method for continuous estimation of clear-sky downwelling longwave
        radiative flux developed using ARM surface measurements.
        Journal of Geophysical Research: Atmospheres, 113, D18206.
        https://doi.org/10.1029/2008JD009936"""
         
        lw_down = self.dataset['downwelling_longwave_irradiation']
        air_temperature = self.dataset['temperature']

        if lw_down.attrs.get("units") != "W/m^2":
            raise ValueError(f"lw_down must have units='W/m^2', has units={lw_down.attrs.get('units')}.")

        if air_temperature.attrs.get("units") != "K":
            raise ValueError(f"air_temperature must have units='K', has units={air_temperature.attrs.get('units')}.")

        valid = (lw_down > 0) & (air_temperature > 0)

        sky_brightness_temperature = (
            lw_down.where(valid) / scipy.constants.Stefan_Boltzmann
        ) ** 0.25

        sky_brightness_temperature = sky_brightness_temperature.rename(
            "sky_brightness_temperature"
        )
        sky_brightness_temperature.attrs = {
            "standard_name": "brightness_temperature",
            "long_name": "effective sky brightness temperature",
            "units": "K",
            "comment": (
                "Equivalent blackbody temperature derived from surface "
                "downwelling longwave irradiance."
            ),
        }

        #todo: decide if i still need this.
        if 0: 
            air_minus_sky_temperature = (
                air_temperature - sky_brightness_temperature
            ).rename("air_minus_sky_temperature")

            air_minus_sky_temperature.attrs = {
                "long_name": (
                    "difference between surface air temperature and "
                    "effective sky brightness temperature"
                ),
                "units": "K",
                "comment": "Calculated as air_temperature - sky_brightness_temperature.",
            }

        return sky_brightness_temperature

    @property
    def apparent_atmospheric_emissivity(self):
        if 'apparent_atmospheric_emissivity' not in self.dataset:
            self.dataset['apparent_atmospheric_emissivity'] = self._calculate_apparent_atmospheric_emissivity()
        return self.dataset['apparent_atmospheric_emissivity']
    
    def _calculate_apparent_atmospheric_emissivity(
            self,
        ) -> xr.DataArray:
        """Calculate apparent broadband atmospheric longwave emissivity.
        
        References
        ----------
        Marty, C., and R. Philipona (2000).
        The clear-sky index to separate clear-sky from cloudy-sky situations
        in climate research.
        Geophysical Research Letters, 27, 2649–2652.
        https://doi.org/10.1029/2000GL011743"""

        lw_down = self.dataset['downwelling_longwave_irradiation']
        air_temperature = self.dataset['temperature']

        if lw_down.attrs.get("units") != "W/m^2":
            raise ValueError(f"lw_down must have units='W/m^2', has units={lw_down.attrs.get('units')}.")

        if air_temperature.attrs.get("units") != "K":
            raise ValueError(f"air_temperature must have units='K', has units={air_temperature.attrs.get('units')}.")

        valid = (lw_down > 0) & (air_temperature > 0)

        apparent_atmospheric_emissivity = (
            lw_down.where(valid)
            / (scipy.constants.Stefan_Boltzmann * air_temperature.where(valid) ** 4)
        ).rename("apparent_atmospheric_emissivity")

        apparent_atmospheric_emissivity.attrs = {
            "long_name": "apparent broadband atmospheric longwave emissivity",
            "units": "1",
            "comment": (
                "Calculated from surface downwelling longwave irradiance "
                "and surface air temperature."
            ),
        }

        return apparent_atmospheric_emissivity