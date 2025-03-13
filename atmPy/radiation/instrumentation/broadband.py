#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 12 12:35:40 2025

@author: hagen
"""
import numpy as np
import re
import scipy as sp
import pandas as pd
import typing



instrument_info = [{'model' : 'SR20-T2', 
                    'inst_type': 'Pyranometer',
                    'manufacturer': 'Hukseflux',},
                   {'model': 'PIR',
                    'inst_type': 'Pyrgeometer',
                    'manufacturer': 'Eppley',},
                   {'model': 'SPN1',
                    'inst_type': 'SPN1',
                    'manufacturer': 'Delta-T Devices'}
                 ]


def parse_calibration_file(filename):
    """
    Probably a very specialized file format ... or is this a common file format for those instruments?
    """
    def find_instrument(searchterm):
        st = searchterm
        res = [a for a in instrument_info if a['model'] == st or a['inst_type'] == st]
        assert(len(res) == 1), f'Found {len(res)} number of instruments for searchterm {st}. Should have found exactly 1'
        return res[0]

    
    with open(filename, 'r', encoding='utf-8') as file:
        txt = file.read()
    
    instruments = re.split(r'\n\s*\n', txt)
       
    cal_list = []
    for inst in instruments:
        lines = inst.strip().split('\n')
        calvalues = {}
        
        line = lines.pop(0)
        if 'Downwelling' in line:
            orientation = 'Downwelling'
        elif 'Upwelling' in line:
            orientation = 'Upwelling'
        else:
            orientation = 'Downwelling' 
        instinfo = find_instrument(line.split()[-1])
    
        
        for line in lines:
            if "SN" in line:
                sn = line.split()[-1]
                calvalues['sn'] = sn
            elif 'Date' in line:
                caldate = line.split(':')[-1].strip()
                calvalues['date_of_calibration'] = caldate
            elif '=' in line:
                try:
                    s, v = line.split('=')
                except ValueError:
                    continue
                try:
                    calvalues[s.strip()] = float(v.strip())
                except ValueError:
                    continue
    
        calvalues['orientation'] = orientation
        calvalues.update(instinfo)
        cal_list.append(calvalues)
    return cal_list


def thermistor_resistance2temperature(r, 
                                      a = 1.0295e-3,
                                      b = 2.391e-4,
                                      c = 1.568e-7,
                                     ):
    r"""
    This function should probably go in a more general corner of the atmPy package ... it might alread exist?
    
    \frac{1}{T} = A + B \ln R + C (\ln R)^3

    Parameters
    ----------
    r: float
        Thermistor resitance in ohm.
    a,b,c: float
        Steinhart-Hart coefficients for a typical 10kΩ thermistor

    Returns
    -------
    Temperatur in Kelvin
    """
    bterm = b*np.log(r)
    cterm = c * np.log(r)**3
    t = 1 / (a + bterm + cterm)
    return t 


class Pyronometer:
    def __init__(self,
                 calibration: typing.Optional[dict] = None, #{'R': 16.392, 'a': -8.8021e-06,'b': 0.00030832, 'c': 0.9974},
                 # company = None,
                 # model = None,
                ):
        """
        Parameters
        ----------
        calibration: dict
            Dictionary with keys: R, a, b, c
            example:
                {'R': 16.392, 'a': -8.8021e-06,'b': 0.00030832, 'c': 0.9974}
        """
        self.calibration = calibration

    def obs2irradiance(self, U, T):
        """
        Calculates irradiance from thermopile voltage and body temperature.

        Parameters
        ----------
        U: float
            Thermopile voltage in mV
        T: float
            Body temperatur in K

        Return
        ------
        irradiance in W/m^2
        """
        
        cal = self.calibration
        T = T - 273.15
        I = U*1e3 / (cal['R'] * ((cal['a'] * T**2) + (cal['b'] * T) + cal['c']))
        
        return I
        



class Pyrgeometer:
    def __init__(self,
                 calibration: dict):
        """
        Parameters
        ----------
        calibration: dict
            Dictionary of calibration coefficients with keys: K1, K2, K3, Kr
            example:
                {'K1': 0.2397, 'K2': 1.0011, 'K3': -4.02, 'Kr': 0.0007044}
        """
        self.calibration = calibration

    def obs2irradiance(self, U, T_dome, T_case):
        """
        Takes the raw observations and derives the irradiance. For details see:
        
        Andreas A, M Dooraghi, A Habte, M Kutchenreiter, I Reda, and M Sangupta. 2018. Solar
        Infrared Radiation Station (SIRS), Sky Radiation (SKYRAD), Ground Radiation (GNDRAD), and
        Broadband Radiometer Station (BRS) Instrument Handbook. Ed. by Robert Stafford, ARM
        Climate Research Facility. DOE/SC-ARM-TR-025. 10.2172/1432706.
        https://www.arm.gov/publications/tech_reports/handbooks/sirs_handbook.pdf

        Parameters
        ----------
        U : float
            Voltage in mV measured by the thermopile. Note Andreas et al. call 
            for micro volts. Since it is more common to measure mV we are considering mV here.
        T_dome : float
            Dome temperature in K.
        T_case : TYPE
            Case temperature in K. Used to estimate the more relevant sensor temperature.

        Returns
        -------
        out : TYPE
            DESCRIPTION.

        """
        U = U * 1e3 
        sig = sp.constants.sigma
        cal = self.calibration
        k4 = cal['Kr']#0.0007044 #K/µV
        T_sensor = T_case + (k4 * U) #T_sensor = T_R
        IR1 = cal['K1'] * U
        IR2 = cal['K2'] * sig * T_sensor**4
        IR3 = cal['K3'] * sig * (T_dome**4 - T_sensor**4)
        IR = IR1 + IR2 + IR3
        out = dict(IR = IR,
                   T_sensor = T_sensor)
        return out
        
        
        
# def read_file(p2f):
#     with open(p2f) as rein:
#         # data = rein.readlines()
#         lines = []
#         for i in range(4):
#             lines.append(rein.readline())
    
#     header = lines[0]
    
#     labels = lines[1]
#     labels = [l.strip().strip('"') for l in labels.split(',')]
    
#     units = lines[2]
#     units = [u.strip().strip('"') for u in units.split(',')]
    
#     smp = lines[3] # no idea what that is and if wee will need it
    
#     dtypes = { a : float for a in labels}
#     dtypes['TIMESTAMP'] = str #pd.Timestamp
    
#     df = pd.read_csv(p2f, skiprows=4, names = labels, index_col=0, 
#                 low_memory=False,
#                 # dtype=dtypes
#                )
    
#     df.index = pd.to_datetime(df.index)
    
#     df.index.name  = 'time'
    
#     # for some reason pandas can figure out all types savely, but the following goes through savely ... makes no sense
#     for dt in df:
#         if df[dt].dtype.name == 'object':
#             df[dt] = df[dt].astype(float)
            
    
#     ds = df.to_xarray()
    
#     # save units in attributs
#     for l,u in zip(labels[1:], units[1:]):
#         ds[l].attrs['unit'] = u
#     ds.attrs['header'] = header
#     return ds