"""This is a collection of classes the are uses as the basis of many retrievals."""

import numpy as np
# import sklearn
import xarray as xr
import pandas as pd
# import atmPy.radiation.retrievals.clearsky as atmcsk
# import matplotlib.pyplot as _plt
from atmPy.opt_imports import matplotlib as mpl
import atmPy.general.measurement_site as atmgms
# import warnings
from atmPy.radiation.retrievals import tiltcorrection
# import pathlib as pl

SOLAR_CONSTANT = 1361.0



class _DatasetRef(object):
    def __init__(self, dataset):
        """A simple class to hold a reference to a dataset. This is used to allow multiple classes to share the same dataset without having to pass it around explicitly."""
        self.dataset = dataset


class SolarIrradiation(object):
    def __init__(self, dataset, site = None, verbose = False, _dataset_ref = None):
        """Base class for solar irradiation datasets. Provides shared functionality and properties for global horizontal, diffuse horizontal, and direct normal irradiation.
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
    def temporal_resolution_minutes(self):
        """Returns the temporal resolution of the dataset in minutes. This is computed as the median difference between consecutive time points."""
        return float((self.dataset.datetime - self.dataset.datetime.shift(datetime = 1)).dropna('datetime').median()/pd.to_timedelta(1, 'm'))

    @property
    def site(self):
        if isinstance(self._site, type(None)):
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
        """Drop variables and update the shared dataset reference."""
        self.dataset = self.dataset.drop_vars(names, errors = errors)
        return self.dataset

    @property
    def mu0(self):
        if 'mu0' not in self.dataset:
            self.dataset['mu0'] = np.cos(self.sun_position.solar_zenith) 
        return self.dataset['mu0']
    
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



class GlobalHorizontalIrradiation(SolarIrradiation):
    def __init__(self, dataset, **kwargs):
        super().__init__(dataset, **kwargs)
        assert('global_horizontal' in dataset), f'global_horizontal variable is missing.'



class DiffuseHorizontalIrradiation(SolarIrradiation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert('diffuse_horizontal' in self.dataset), f'diffuse_horizontal variable is missing.'


class DirectNormalIrradiation(SolarIrradiation):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert(('direct_normal' in self.dataset) or ('direct_horizontal' in self.dataset)), f'direct_normal or alternatively direct_horizontal variable is missing.'

        if 'direct_horizontal' in self.dataset and 'direct_normal' not in self.dataset:
            self.direct_normal_from_direct_horizontal()


    def direct_normal_from_direct_horizontal(self):
        """Converts direct horizontal irradiance to direct normal irradiance. 'direct_normal' will be added to the dataset if it is not already present. """
        if 'direct_normal' not in self.dataset:
            ds = self.dataset
            ds['direct_normal'] = ds.direct_horizontal / xr.DataArray(np.cos(self.sun_position.solar_zenith))
            ds.direct_normal.attrs.update(ds.direct_horizontal.attrs)
            ds.direct_normal.attrs["long_name"] += ' (calculated from direct_horizontal and solar zenith angle)'
        return ds['direct_normal']
    

class CombinedGlobalDiffuseDirect(SolarIrradiation):
    def __init__(self, dataset, **kwargs):
        """
        A class that combines the three irradiation components into one dataset and provides some specific functions
        Parameters
        ----------
        dataset : xr.Dataset
            A dataset that contains the three components as variables named 'global_horizontal', 'diffuse_horizontal', 'direct_normal'
            and the coordinates 'datetime' and 'channel'
        -----------
        2024-06-10, HTelg   
        """

        super().__init__(dataset, **kwargs)
        self._global_horizontal_irradiation = None
        self._diffuse_horizontal_irradiation = None
        self._direct_normal_irradiation = None
        self._kwargs = kwargs

    
    # @property
    # def clearsky_global_irradiation_powerlow_fit_params(self):
    #     """Performes a empirical power law fit to the clear sky global irradiance
    #     and returns the fit parameters."""

    #     if 'clearsky_global_irradiation_powerlow_fit_params' not in self.dataset:
    #         if self.verbose:
    #             print('performing powerlow fit for clear sky global irradiance')
    #         params = atmcsk.fit_global_powerlaw_mu0(self.mu0, self.dataset.global_horizontal, self.mask_clear_sky_radflux, 
    #                                                     mu0_min = self.get_attr('mu0_min'),
    #                                                     min_points = 100)
    #         if isinstance(params, type(None)) and self.verbose:
    #             print('fit failed, probably not enough clear sky points')
                
    #         self.dataset['clearsky_global_irradiation_powerlow_fit_params'] = params
    #     return self.dataset['clearsky_global_irradiation_powerlow_fit_params']

    # @property
    # def clearsky_diffuse_irradiation_powerlow_fit_params(self):
    #     """Performes a empirical power law fit to the clear sky diffuse irradiance
    #     and returns the fit parameters."""

    #     if 'clearsky_diffuse_irradiation_powerlow_fit_params' not in self.dataset:
    #         if self.verbose:
    #             print('performing powerlow fit for clear sky diffuse irradiance')
    #         dt_in_m = np.median(self.dataset.datetime.values[1:]-self.dataset.datetime.values[:-1])/pd.to_timedelta(1,'m')

    #         params = atmcsk.fit_global_powerlaw_mu0(self.mu0, self.dataset.diffuse_horizontal, self.mask_clear_sky_radflux, 
    #                                                     mu0_min = self.get_attr('mu0_min'),
    #                                                     min_points = dt_in_m * 100)
    #         if isinstance(params, type(None)) and self.verbose:
    #             print('fit failed, probably not enough clear sky points')
                
    #         self.dataset['clearsky_diffuse_irradiation_powerlow_fit_params'] = params
    #     return self.dataset['clearsky_diffuse_irradiation_powerlow_fit_params']
    


    # def plot_clearsky_global_irradiation_powerlow_fit(self, ax = None):
    #     params = self.clearsky_global_irradiation_powerlow_fit_params
    #     if any(params.isnull()):
    #         print('Not enough clearsky points for valid fit.')
    #         return None, None
    #     if isinstance(ax, type(None)):
    #         f, a= mpl.pyplot.subplots()    
    #     else:
    #         a = ax
    #         f = a.get_figure()
        
    #     A = params.to_pandas().a
    #     b = params.to_pandas().b
    #     fit = A * np.power(self.mu0, b)
    #     fit.plot(ax = a, label = 'fit')
    #     return f,a


    @property
    def global_horizontal_irradiation(self):
        if isinstance(self._global_horizontal_irradiation, type(None)):
            self._global_horizontal_irradiation = GlobalHorizontalIrradiation(
                self.dataset,
                _dataset_ref = self._dataset_ref,
                **self._kwargs,
            )
        return self._global_horizontal_irradiation
    
    @property
    def diffuse_horizontal_irradiation(self):
        if isinstance(self._diffuse_horizontal_irradiation, type(None)):
            self._diffuse_horizontal_irradiation = DiffuseHorizontalIrradiation(
                self.dataset,
                _dataset_ref = self._dataset_ref,
                **self._kwargs,
            )
        return self._diffuse_horizontal_irradiation

    @property
    def direct_normal_irradiation(self):
        if isinstance(self._direct_normal_irradiation, type(None)):
            self._direct_normal_irradiation = DirectNormalIrradiation(
                self.dataset,
                _dataset_ref = self._dataset_ref,
                **self._kwargs,
            )
        return self._direct_normal_irradiation
    

    def plot_overview(self, 
                    #   channel = 500, 
                      ax = None,
                      show_clearsky = False, 
                      apply_mask_clear_sky = False,
                      show_sunelevation = False,
                      plot_kwargs = {},):
        
        if isinstance(ax, type(None)):
            f, a= mpl.pyplot.subplots()    
        else:
            a = ax
            f = a.get_figure()
        
        self.global_horizontal_irradiation
        self.diffuse_horizontal_irradiation
        self.direct_normal_irradiation
        
        dssel = self.dataset

        dssel.global_horizontal.plot(ax = a, label = 'global_horizontal')
        if apply_mask_clear_sky:
            g = a.get_lines()[-1]
            col = g.get_color()
            g.set_alpha(0.5)
            g.set_label('_nolegend_')
            dssel.global_horizontal.where(self.mask_clear_sky_radflux).plot(ax = a, label = 'global_horizontal', color = col)
        dssel.diffuse_horizontal.plot(ax = a, label = 'diffuse_horizontal')
        if apply_mask_clear_sky:
            g = a.get_lines()[-1]
            col = g.get_color()
            g.set_alpha(0.5)
            g.set_label('_nolegend_')
            dssel.diffuse_horizontal.where(self.mask_clear_sky_radflux).plot(ax = a, label = 'diffuse_horizontal', color = col)
        dssel.direct_normal.plot(ax = a, label = 'direct')
        if apply_mask_clear_sky:
            g = a.get_lines()[-1]
            col = g.get_color()
            g.set_alpha(0.5)
            g.set_label('_nolegend_')
            dssel.direct_normal.where(self.mask_clear_sky_radflux).plot(ax = a, label = 'direct_normal', color = col)
        
        if show_clearsky:
            self.clearsky_global_horizontal.plot(ax = a, ls = '--', label = 'clearsky_global_horizontal')
            self.clearsky_diffuse_horizontal.plot(ax = a, ls = '--', label = 'clearsky_diffuse_horizontal')

        if show_sunelevation:
            at = a.twinx()
            np.rad2deg(self.sun_position.elevation).plot(ax = at, color = 'black', ls = '--')
            # at.set_ylim(top = 0.9, bottom = 0)
        
        # a.set_xlim(left = pd.to_datetime('20220103 14:00:00'))
        a.grid()
        a.legend()
        return f,a




    def apply_tilt_correction(self, sensor_response_time=None):
        # call irradiance properties to make sure they are initialized and data is present in the dataset
        # self.direct_normal_irradiation  
        self.sun_position
        ds = tiltcorrection.apply_tilt_correction(self.dataset, sensor_response_time=sensor_response_time)
        out = CombinedGlobalDiffuseDirect(ds)
        return out

    apply_tilt_correction.__doc__ = tiltcorrection.apply_tilt_correction.__doc__

    def convert2RadFlux(self):
        """Convert the dataset to a atmPy.radiation.radflux.lab.RadFlux object. """
        from atmPy.radiation.radflux.lab import RadFlux

        return RadFlux(self.dataset, verbose = self.verbose, site = self.site)
