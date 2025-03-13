from . import broadband

import xarray as xr

class RadSys:
    def __init__(self, calibration: list, 
                 # sn_sw_up: typing.Optional[str] = None,
                 # sn_sw_down: typing.Optional[str] = None,
                 # sn_lw_up: typing.Optional[str] = None,
                 # sn_lw_down: typing.Optional[str] = None,
                 # sn_spn1: typing.Optional[str] = None,
                 verbose: bool = True, 
                ):
        """
        Example calibration variable
        ----------------------------
        [{'sn': '15337',
          'date_of_calibration': '4/29/24',
          'R': 17.033,
          'a': -8.8146e-06,
          'b': 0.00035343,
          'c': 0.9965,
          'orientation': 'Downwelling',
          'model': 'SR20-T2',
          'inst_type': 'Pyranometer',
          'manufacturer': 'Hukseflux'},
         {'sn': '15340',
          'date_of_calibration': '4/29/24',
          'R': 16.392,
          'a': -8.8021e-06,
          'b': 0.00030832,
          'c': 0.9974,
          'orientation': 'Upwelling',
          'model': 'SR20-T2',
          'inst_type': 'Pyranometer',
          'manufacturer': 'Hukseflux'},
         {'sn': '37336F3',
          'date_of_calibration': '5/21/24',
          'K1': 0.2397,
          'K2': 1.0011,
          'K3': -4.02,
          'Kr': 0.0007044,
          'orientation': 'Downwelling',
          'model': 'PIR',
          'inst_type': 'Pyrgeometer',
          'manufacturer': 'Eppley'},
         {'sn': '37338F3',
          'date_of_calibration': '5/21/24',
          'K1': 0.23904,
          'K2': 1.004,
          'K3': -4.7,
          'Kr': 0.0007044,
          'orientation': 'Upwelling',
          'model': 'PIR',
          'inst_type': 'Pyrgeometer',
          'manufacturer': 'Eppley'},
         {'sn': 'A4127',
          'C': 1.0,
          'orientation': 'Downwelling',
          'model': 'SPN1',
          'inst_type': 'SPN1',
          'manufacturer': 'Delta-T Devices'}]
        """
        self.verbose = verbose
        self.calibration = calibration.copy()

        def _find_in_cal(orientation, inst_type):
            cal = [cal for cal in self.calibration if (cal["orientation"] == orientation) and (cal["inst_type"] == inst_type)]
            if len(cal) == 0:
                raise ValueError(f'No calibration found with orientation:{orientation} and inst_type:{inst_type}.')
            elif len(cal) > 1:
                raise ValueError(f'More than 1 calibration found with orientation:{orientation} and inst_type:{inst_type}.')
            cal = cal[0]
            if self.verbose:
                print(f'found cal:\n{cal}')
            return cal
            

        cal = _find_in_cal('Downwelling', 'Pyranometer')
        self.sw_down = broadband.Pyronometer(calibration = cal)
        
        cal = _find_in_cal('Upwelling', 'Pyranometer')
        self.sw_up = broadband.Pyronometer(calibration = cal)

        cal = _find_in_cal('Downwelling', 'Pyrgeometer')
        self.lw_down = broadband.Pyrgeometer(cal)

        cal = _find_in_cal('Upwelling', 'Pyrgeometer')
        self.lw_up = broadband.Pyrgeometer(cal)

        # cal = _find_in_cal('Downwelling', 'Pyranometer')
        # self.sw_down = Pyronometer(cal)

        # cal = _find_in_cal('Downwelling', 'Pyranometer')
        # self.sw_down = Pyronometer(cal)

        # cal = _find_in_cal('Downwelling', 'Pyranometer')
        # self.sw_down = Pyronometer(cal)



        # if (sn := sn_sw_up) is not None:
        #     if not isinstance(sn, str):
        #         raise TypeError(f"kwarg must be a string, but got {type(sn).__name__}")
        #     cal = [cal for cal in calibration if cal[sn] == sn] 
        #     if len(cal) == 0:
        #         raise ValueError(f'Serial no {sn} not found in calibration.')
        #     elif len(cal) > 1:
        #         raise ValueError(f'Serial no {sn} appears more than once in calibration.')
                
        #     cal = cal[0]
        #     self.sw_up = Pyronometer(calibration = cal)

    def obs2radiation(self,observation):
        ds_obs = observation
        ds_out = xr.Dataset()
        
        ds_out['inst_down_short_hemisp_case_temp'] = broadband.thermistor_resistance2temperature(ds_obs.inst_down_short_global_case_resist * 1e3)
        ds_out['inst_down_short_hemisp'] = self.sw_down.obs2irradiance(ds_obs.inst_down_short_tp, ds_out['inst_down_short_hemisp_case_temp'])
        
        ds_out['inst_up_short_hemisp_case_temp'] = broadband.thermistor_resistance2temperature(ds_obs.inst_up_short_case_resist * 1e3)
        ds_out['inst_up_short_hemisp'] = self.sw_up.obs2irradiance(ds_obs.inst_up_short_tp, ds_out['inst_up_short_hemisp_case_temp'])
        
        t_case = broadband.thermistor_resistance2temperature(ds_obs.inst_up_long_hemisp_case_resist * 1e3)
        ds_out['inst_up_long_hemisp_case_temp'] = t_case
        t_dome = broadband.thermistor_resistance2temperature(ds_obs.inst_up_long_hemisp_dome_resist * 1e3)
        ds_out['inst_up_long_hemisp_dome_temp'] = t_dome
        
        out = self.lw_up.obs2irradiance( ds_obs.inst_up_long_hemisp_tp, t_dome, t_case)
        ds_out['inst_up_long_hemisp'] = out['IR']
        
        t_case = broadband.thermistor_resistance2temperature(ds_obs.inst_down_long_hemisp_case_resist * 1e3)
        ds_out['inst_down_long_hemisp_case_temp'] = t_case
        t_dome = broadband.thermistor_resistance2temperature(ds_obs.inst_down_long_hemisp_dome_resist * 1e3)
        ds_out['inst_down_long_hemisp_dome_temp'] = t_dome
        
        out = self.lw_down.obs2irradiance(ds_obs.inst_down_long_hemisp_tp, t_dome, t_case)
        ds_out['inst_down_long_hemisp'] = out['IR']
        
        return ds_out

        
            
                