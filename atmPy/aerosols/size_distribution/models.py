#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 12:47:56 2022

@author: hagen
"""

import pandas as pd
import atmPy.aerosols.size_distribution.sizedistribution as sd
import matplotlib.pyplot as plt


class Model(object):
    def __init__(self, model_params, name = None):
        
        self._model_params = model_params
        self._diameter_range = [2, 1e5] # nanometer
        self._number_of_diameters = 100
        
        self._reset()
        
    def _reset(self):
        self._size_distribution = None
    
    @property
    def model_params(self):
        return self._model_params
    
    @model_params.setter
    def model_params(self,value):
        self._model_params = value
        self._reset()
        return
    
    @property
    def diameter_range(self):
        return self._diameter_range
    
    @diameter_range.setter
    def diameter_range(self,value):
        self._diameter_range = value
        self._reset()
        return
    
    @property
    def number_of_diameters(self):
        return self._number_of_diameters
    
    @number_of_diameters.setter
    def number_of_diameters(self,value):
        self._number_of_diameters = value
        self._reset()
        return
    
    
    @property
    def size_distribution(self):
        if isinstance(self._size_distribution,type(None)):
            dists = []
            for idx, row in self.model_params.iterrows():
                distt = sd.simulate_sizedistribution(diameter=self.diameter_range,
                                            numberOfDiameters=self.number_of_diameters,
                                            centerOfAerosolMode=row.diameter,
                                            widthOfAerosolMode=row.sig,
                                            numberOfParticsInMode=row.number)
                distt = distt.convert2dNdDp()
                dists.append(distt)
            #     break
            
            dist = dists[0].copy()
            for distt in dists[1:]:
                dist += distt
            dist = dist.convert2dNdlogDp()
            self._size_distribution = dist
        return self._size_distribution
  

class SeignfeldAndPandas(object):
    def __init__(self):
        """
        From Siegnfeld and Pandas, third edition, Page 343. Table 8.3.
        Note, not all of the size distributions that are depected on the 
        following pages can be reproduced. However, since I was able to 
        reproduce some, eg. rural and free-troposphere I assume my code is 
        correct. In addition to discrepancies in the general shape there is a 
        slight deviation in total number too (even for rural and free 
        troposphere). Since the number concentration function produces the
        right number I again have confident in my code. (@myself: tests can be 
        found here: projects/16_closure_of_arm_data/uncertainties/aerosol_models/continental/uncertainties_scattering_only.ipynb)

        Returns
        -------
        None.

        """
        self._model_list = self._get_model_list()
        for mod in self._model_list:
            model = Model(mod['model_params'], name = mod['zone'])
            setattr(self, mod['zone'], model)
            mod['model'] = model
        
    def plot_sizedistributions(self, ax = None, moment ='number', **kwargs):
        if isinstance(ax, type(None)):
            f,a = plt.subplots()
        else:
            a = ax
        for mod in self._model_list:
            model = mod['model']
            dist = model.size_distribution
            if moment == 'number':
                dist = dist.convert2dNdlogDp()
            elif moment == 'surface':
                dist = dist.convert2dSdlogDp()
            elif moment == 'volume':
                dist = dist.convert2dVdlogDp()
            dist.plot(ax = a, label = mod['zone'], **kwargs)
        a.set_yscale('log')
        a.set_ylim(1e-3, 2e4)
        a.legend()
        return a
    
    def _get_model_list(self):
        models = [] # a collection of data for all models

        ## Urban
        
        mode1 = {'number': 7100 ,'diameter': 0.0117,   'sig': 0.232}
        mode2 = {'number': 6320,'diameter': 0.0373,   'sig': 0.25}
        mode3 = {'number': 960,'diameter': 0.151,   'sig': 0.204}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        models.append({'zone': 'urban', 'model_params': model_params})
        
        ## Marine
        
        mode1 = {'number':  133,'diameter': 0.008,   'sig': 0.657}
        mode2 = {'number': 66.6,'diameter': 0.266,   'sig': 0.21}
        mode3 = {'number': 3.1,'diameter': 0.58,   'sig': 0.396}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        
        models.append({'zone': 'marine', 'model_params': model_params})
        
        ## Rural
        
        mode1 = {'number':  6650,'diameter': 0.015,   'sig': 0.225}
        mode2 = {'number': 147,'diameter': 0.054,   'sig': 0.557}
        mode3 = {'number': 1990,'diameter': 0.084,   'sig': 0.266}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        
        models.append({'zone': 'rural', 'model_params': model_params})
        
        ## remote continental
        
        mode1 = {'number':  3200,'diameter': 0.02,   'sig': 0.161}
        mode2 = {'number': 2900,'diameter': 0.116,   'sig': 0.217}
        mode3 = {'number': 0.3,'diameter': 1.8,   'sig': 0.380}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        
        models.append({'zone': 'continental', 'model_params': model_params})
        
        ## free troposphere
        
        mode1 = {'number':  129,'diameter': 0.007,   'sig': 0.645}
        mode2 = {'number': 59.7,'diameter': 0.25,   'sig': 0.253}
        mode3 = {'number': 63.5,'diameter': 0.52,   'sig': 0.425}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        
        models.append({'zone': 'free_trop', 'model_params': model_params})
        
        ## Polar
        
        mode1 = {'number':  21.7,'diameter': 0.138,   'sig': 0.245}
        mode2 = {'number': 0.186,'diameter': 0.75,   'sig': 0.3}
        mode3 = {'number': 3e-4,'diameter': 8.6,   'sig': 0.291}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        
        models.append({'zone': 'polar', 'model_params': model_params})
        
        ## desert
        
        mode1 = {'number':  726,'diameter': 0.002,   'sig': 0.247}
        mode2 = {'number': 114,'diameter': 0.038,   'sig': 0.770}
        mode3 = {'number': 0.178,'diameter': 21.6,   'sig': 0.438}
        model_params = pd.DataFrame([mode1, 
                                      mode2, 
                                      mode3])
        model_params.loc[:,'diameter'] *= 1e3
        
        models.append({'zone': 'desert', 'model_params': model_params})
        return models
        