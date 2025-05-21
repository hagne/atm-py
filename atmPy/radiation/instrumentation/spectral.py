# -*- coding: utf-8 -*-
import xarray as xr
import numpy as np
import pathlib as pl
import atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab as atmcal
import atmPy.radiation.retrievals.spectral_irradiance as atmsi
from .. import solar as atmsol


class Mfr():
    def __init__(self,
                 spectral_calibration = None,
                 logger_calibration = None,
                 head_calibration = None,
                 cosine_responds = None):
        
        self.spectral_calibration = spectral_calibration
        self.logger_calibration = logger_calibration
        if np.array_equal(logger_calibration, head_calibration):
            print('logger equal head calibration')
            self.head_calibration = self.logger_calibration.copy()
        else:
            self.head_calibration = head_calibration
        self.cosine_responds = cosine_responds
        
        self._data_class = atmsi.GlobalHorizontalIrradiation
        
        self._cosine_calibration_diffuse = None
        
    @property
    def spectral_calibration(self):
        return self._spectral_calibration
    
    @spectral_calibration.setter 
    def spectral_calibration(self, value):
        if isinstance(value, xr.Dataset):
            self._spectral_calibration = value
        elif isinstance(value, (str, pl.Path)):
            self._spectral_calibration = atmcal.read_mfrsr_cal(value)
        else:
            raise TypeError(f'Unknown type for spectral_calibration: {type(value)}')
        assert('statistics' in self._spectral_calibration.variables), "I don't think the calibration file is a spectral calibration file for an MFR(SR). 'statistics' varible is missing"
        return 
    
    @property
    def logger_calibration(self):
        return self._logger_calibration
    
    @logger_calibration.setter
    def logger_calibration(self,value):
        if isinstance(value, xr.Dataset):
            self._logger_calibration = value
        elif isinstance(value, (str, pl.Path)):
            self._logger_calibration = atmcal.read_factory_cal(value)
        else:
            raise TypeError(f'Unknown type for logger_calibration: {type(value)}')
            
        return
    
    @property
    def head_calibration(self):
        return self._head_calibration
    
    @head_calibration.setter 
    def head_calibration(self, value):
        if isinstance(value, xr.Dataset):
            self._head_calibration = value
        elif value is None:
            self._head_calibration = None
        elif isinstance(value, (str, pl.Path)):
            if pl.Path(value).suffix == '.nc':
                self._head_calibration = xr.open_dataset(value)
            else:
                self._head_calibration = atmcal.read_mfrsr_cal_responsivity(value)
        else:
            raise TypeError(f'Unknown type for head_calibration: {type(value)}')
            
        return
    
    @property
    def cosine_responds(self):
        return self._cosine_responds
    
    @cosine_responds.setter 
    def cosine_responds(self, value):
        if isinstance(value, xr.Dataset):
            self._cosine_responds = value
        elif isinstance(value, (str, pl.Path)):
            try:
                self._cosine_responds = atmcal.read_mfrsr_cal_cos(value)
            except AttributeError:
                try:
                    self._cosine_responds = atmcal.read_mfrsr_cal_cos(value, 
                                                                      broadband_col_name='Termopile', 
                                           )
                except AttributeError:
                    try:
                        self._cosine_responds = atmcal.read_mfrsr_cal_cos(value, 
                                                                          broadband_col_name='300')
                    except:
                        raise
        elif value is None:
            self._cosine_responds = None
        else:
            raise TypeError(f'Unknown type for cosine_responds: {type(value)}')
            
        return 
    
    def _apply_calibration_spectral(self, si):
        """
        This will assign the nominal channel wavelength and provide the exact
        channel central wavelengths.

        Parameters
        ----------
        si : SpectralIrradiance class or subclass
            Use atmPy.data_archives.NOAA_ESRL_GMD_GRAD.cal_facility.lab.read_mfrsr_cal 
            to open calibration file and use return as input here.

        Returns
        -------
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.

        """
        
        calibration = self.spectral_calibration
        
        ds = si.dataset.copy()     
        
        if ds.channel.shape[0] == 7:            
            ds = ds.isel(channel = slice(1,None)) # this gets rid of the broadband channel
        if 'global_horizontal' not in ds.variables:
            ds = ds.rename_vars({'alltime': 'global_horizontal'})
        
        
        ds['channel'] = calibration.channel
        ds['channel_wavelength'] = calibration.statistics.sel(stats = 'CENT', drop=True)
        ds.attrs['calibrated_spectral'] = 'True'
        ds = ds.sortby('channel')
        return self._data_class(ds) #returns the same class, allows for application to all subclasses
    
    def _apply_calibration_responsivity(self, si, calibration, 
                                       # varname_responsivity_spectral = 'responsivity_spectral',
                                       # varname_dark_signal_spectral = 'dark_signal_spectral',
                                       # ignore_has_been_applied_error = False
                                       ):
        """
        This will calibrate for amplifier responsivity. This is sometimes
        applied multiple times, e.g. in MFR-type instruments where we have head
        and datalogger sensitivity. If executed multiple time the "has been 
        applied error" will trigger. Make sure to set ignore_has_been_applied_error
        to True

        Parameters
        ----------
        calibration : TYPE
            DESCRIPTION.
         : TYPE
            DESCRIPTION.

        Returns
        -------
        TYPE
            DESCRIPTION.
        TYPE
            DESCRIPTION.

        """

        
        # if not ignore_has_been_applied_error:
        #     if 'calibration_dark_signal' in self.dataset.attrs:
        #         assert(self.dataset.attrs['calibration_dark_signal'] != 'True'), 'Responds calibration already applied'        
        #     if 'calibration_responds' in self.dataset.attrs:
        #         assert(self.dataset.attrs['calibration_responds'] != 'True'), 'Responds calibration already applied'

        assert(isinstance(calibration, xr.Dataset))
        assert('dark_signal_spectral' in calibration.variables), "I don't think the calibration file is a spectral calibration file for an MFR(SR). 'dark_signal_spectral' varible is missing"
        
        ds = si.dataset.copy() 
        #### global horizontal
        
        #### dark signal
        da = ds.global_horizontal - calibration.dark_signal_spectral 
        #### responsivity
        da = da / calibration.responsivity_spectral
        
        # da.attrs['unit'] = 'W * m^-2 * nm'
        # da.attrs['calibration_responds'] = 'True'
        # da.attrs['calibration_dark_signal'] = 'True'
        
        ds['global_horizontal'] = da
        
        #### diffuse horizontal
        if 'diffuse_horizontal' in ds.variables:
            da = ds.diffuse_horizontal - calibration.dark_signal_spectral
            da = da / calibration.responsivity_spectral
            
            # da.attrs['unit'] = 'W * m^-2 * nm'
            # da.attrs['calibration_dark_signal'] = 'True'
            # da.attrs['calibration_responds'] = 'True'
            
            ds['diffuse_horizontal'] = da
            
        # ds.attrs['calibration_dark_signal'] = 'True'
        # ds.attrs['calibration_responds'] = 'True'
        return self._data_class(ds)
    
    def _apply_calibration_logger(self,si):
        ### logger gain calibration
        calibration = self.logger_calibration.copy()
        if 'logger_gain' in calibration.variables:
            "in case of factory calibratios both head and logger calibrations are in the same file"
            calibration = calibration.rename({'logger_gain':'responsivity_spectral',
                                              'logger_darksignal': 'dark_signal_spectral'})
        else:
            assert(False), 'programming required'
            
        # mfr = mfr.apply_calibration_responsivity(mfr_cal_logger, 
        #                                                            varname_responsivity_spectral = varname_responsivity_spectral,
        #                                                            varname_dark_signal_spectral = varname_dark_signal_spectral,
        #                                                        # ignore_has_been_applied_error = True,
        #                                                       ) #logger calibration from factory cal ==> counts to voltage
        # mfr_logger = mfr
        
        return self._apply_calibration_responsivity(si, calibration)
    
    def _apply_calibration_head(self,si):
        calibration = self.head_calibration
        if calibration is None:
            return si
        if calibration.attrs['product_name'] == 'side_by_side_calibration':
            ds = si.dataset.copy()
            ds['global_horizontal'] = ds.global_horizontal * calibration.relative_MFR_head_responsivity
            return self._data_class(ds)
        else:
            if 'head_gain' in calibration.variables:
                "in case of factory calibratios both head and logger calibrations are in the same file"                    
                calibration = calibration.rename({'head_gain':'responsivity_spectral',
                                                  'head_darksignal': 'dark_signal_spectral'})
            # else:
            #     assert(False), 'programming required'
            return self._apply_calibration_responsivity(si, calibration)
    
    @property
    def cosine_calibration_diffuse(self):
        #### for diffuse or global in case of an MFR
        if self._cosine_calibration_diffuse is None:
            cal_angle = 45
            ew = self.cosine_responds.spectral_EW.interp(Angle = [cal_angle, -cal_angle]).sum(dim = 'Angle') / 2 
            ns = self.cosine_responds.spectral_NS.interp(Angle = [cal_angle, -cal_angle]).sum(dim = 'Angle') / 2 
            cal = (ew + ns) / 2
            self._cosine_calibration_diffuse = 1/cal
        return self._cosine_calibration_diffuse
   
    def _apply_calibration_cosine(self, si):
        # if 'clalibration_cosine' in self.dataset.attrs:
        #     assert(self.dataset.attrs['clalibration_cosine'] != 'True'), 'Responds calibration already applied'  
        if self.cosine_responds is None:
            return si
        ds = si.dataset.copy()
        ds['global_horizontal'] = ds.global_horizontal * self.cosine_calibration_diffuse
        ds['cosine_cal_coeff_diffuse'] = self.cosine_calibration_diffuse
        return self._data_class(ds)
        
    def get_cosine_calibration_direct(self, si):
        """
        This depends on the actual position of the sun, therefore it is not 
        intrinsic to the instrument like the cosine correction for the diffuse,
        but needs the actual data to calculate the sun position.

        Parameters
        ----------
        si : TYPE
            DESCRIPTION.

        Returns
        -------
        None.

        """
        calibration = self.cosine_responds

        #### - direct
        sp = atmsol.SolarPosition(si.sun_position.azimuth, np.pi/2 - si.sun_position.elevation, unit = 'rad')
        
        # NS
        # interpolate the cosine respond with the particular angles resulting from the projetion
        # This results in :
        #     * calibration value as a function of time
        
        da = calibration.spectral_NS.interp(Angle = np.rad2deg(sp.projectionNS_angle) - 90)
        cos_cal_NS = da.rename({'Angle': 'datetime'}).assign_coords(datetime = ('datetime', si.sun_position.index))
        
        
        # The calibration value needs to be normalized with the relevant component of the solar radiation
        # With other words how much light is actually comming this way?
        cos_cal_NS_norm = cos_cal_NS * xr.DataArray(sp.projectionNS_norm)
        
        # Do the same for EW
        
        da = calibration.spectral_EW.interp(Angle = np.rad2deg(sp.projectionEW_angle) - 90)
        cos_cal_EW = da.rename({'Angle': 'datetime'}).assign_coords(datetime = ('datetime', si.sun_position.index))
        cos_cal_EW_norm = cos_cal_EW * xr.DataArray(sp.projectionEW_norm)
        
        # Sum NS and EW
        cos_cal_sum = cos_cal_EW_norm + cos_cal_NS_norm
        
        # Divide by the sum of **norms** (Not the calibration value! As we are dealing with vectors the sum is not automatically num
        sumofnorm = sp.projectionEW_norm + sp.projectionNS_norm
        cos_cal_sum_nom = cos_cal_sum / xr.DataArray(sumofnorm)
        cos_cal_sum_nom = 1/cos_cal_sum_nom
        return cos_cal_sum_nom
        
    
    def raw2calibrated(self, data):
        if isinstance(data, (str, pl.Path)):
            ds = xr.open_dataset(data)
        else:
            raise TypeError(f'Unknown type for data: {type(data)}')
        
        # make sure that direct obs are removed as they are produced on the fly and only result in confusion
        ds = ds.drop_vars(['direct_normal', 'direct_horizontal'], errors = 'ignore')
        
        si = self._data_class(ds)
        self.data_raw = si
        
        si = self._apply_calibration_spectral(si)
        self.data_cal_spec = si
        
        si = self._apply_calibration_logger(si)
        self.data_cal_spec_log = si
        
        si = self._apply_calibration_head(si)
        self.data_cal_spec_log_head = si
        
        si = self._apply_calibration_cosine(si)
        self.data_cal_spec_log_head_cos = si
        
        return si
            
class Mfrsr(Mfr):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._data_class = atmsi.CombinedGlobalDiffuseDirect
        
    def _apply_calibration_cosine(self, si):
        # if 'clalibration_cosine' in self.dataset.attrs:
        #     assert(self.dataset.attrs['clalibration_cosine'] != 'True'), 'Responds calibration already applied'  
        if self.cosine_responds is None:
            return si
        
        ds = si.dataset.copy()
        
        ds['direct_horizontal'] = (ds.global_horizontal - ds.diffuse_horizontal)
        #### - diffuse
        ds['diffuse_horizontal'] = ds.diffuse_horizontal * self.cosine_calibration_diffuse        
        cos_cal_dir = self.get_cosine_calibration_direct(si)
        ds['direct_horizontal'] = ds.direct_horizontal * cos_cal_dir
        ds['cosine_calibraion_direct'] = cos_cal_dir
        # ds['direct_normal'] = ds.direct_horizontal / xr.DataArray(np.sin(si.sun_position.elevation))
        
        ds['solar_zenith_angle'] = 90 - np.rad2deg(si.sun_position.elevation) 
        ds['solar_azimuth_angle'] = np.rad2deg(si.sun_position.azimuth)
        #### - compose global based on cosine corrected direct and diffuse
        ds['global_horizontal'] = ds.direct_horizontal + ds.diffuse_horizontal
            
        ds.attrs['clalibration_cosine'] = 'True'
        out = self._data_class(ds) #returns the same class, allows for application to all subclasses
        # if 'diffuse_horizontal' in ds.variables:
        #     out.tp_cos_cal_sum = cos_cal_sum
        #     out.tp_cos_cal_EW_norm = cos_cal_EW_norm
        #     out.tp_cos_cal_NS_norm = cos_cal_NS_norm
        #     out.tp_cos_cal_sum_nom = cos_cal_sum_nom
        return out
    
    def raw2calibrated(self, data):
        si = super().raw2calibrated(data)
        if 'direct_horizontal' not in si.dataset:
            si.dataset['direct_horizontal'] = (si.dataset.global_horizontal - si.dataset.diffuse_horizontal)
        sialt = self.data_cal_spec_log_head # there is a chance that this instance already has the sunposition calculated, avoids double execution
        si.dataset['direct_normal'] = si.dataset.direct_horizontal / xr.DataArray(np.sin(sialt.sun_position.elevation))
        return si