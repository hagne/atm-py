#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 12 12:47:56 2022

@author: hagen
"""

import pandas as pd
import atmPy.aerosols.size_distribution.sizedistribution as sd
import atmPy.aerosols.size_distribution.diameter_binning as db
import matplotlib.pyplot as plt
import numpy as np



class Model(object):
    def __init__(self, model_params, name = None, diameter_range = [2, 2e5]):
        
        self._model_params = model_params
        self._diameter_range = diameter_range # nanometer
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
    def __init__(self, **kwargs):
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
            model = Model(mod['model_params'], name = mod['zone'], **kwargs)
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
    
class SatellitAerosolModelsABI(object):
    def __init__(self, aod=0.1, diameter_range = [1e2, 2e4, 100]):
        """
        WARNING: 
            I am pretty sure there is a normalization problem, since I 
            substitude um with nm and r with d without doing proper normailzing.
            Also, the model is created for the natural logarithm while atmPy 
            assumes a log_10! This requires a further normalizaion (see 
            Seignfeld & Pandis). These things do not affect general shape but 
            will need to be addressed if absolute values are considered.
            
        Aerosol models used by the ABI aerosol optical depth retrieval. 
        From:
            GOES-R Advanced Baseline Imager (ABI) Algorithm Theoretical Basis
            Document For Suspended Matter/Aerosol Optical Depth and Aerosol 
            Size Parameter
            https://www.goes-r.gov/resources/docs.html
        The aerosol models of satellite retrievals do not strictly follow my
        model class so this is not inheriting Model at this point ... maybe
        later?
        
        Parameters
        ----------
        aod : float
            Aerosol optical depth. The exact aerosol model depends on the aerosol optical depth.
        diameter_range : array-like, optional
            Diameter range, in nanometer, and number of points the model is created for. The default is [1e2, 1e4, 100].

        Returns
        -------
        None.

        """
        models = pd.DataFrame([   {'model': 'generic', 'mode': 'fine',   'rv': 0.145 , 'rv_scale': 0.0203, 'sig': 0.3738, 'sig_scale': 0.1365, 'Cv': .1642 , 'Cv_scale': 0.7747, 'n_r': 1.43, 'n_i': 0.008, 'n_scale': 0.002},
                          {'model': 'generic', 'mode': 'coarse', 'rv': 3.1007, 'rv_scale': 0.3364, 'sig': 0.7292, 'sig_scale': 0.098 , 'Cv': 0.1482, 'Cv_scale': 0.6846, 'n_r': 1.43, 'n_i': 0.008, 'n_scale': 0.002},
                          {'model': 'urban',   'mode': 'fine',   'rv': 0.1604, 'rv_scale': 0.434,  'sig': 0.3642, 'sig_scale': 0.1529, 'Cv': 0.1718, 'Cv_scale': 0.8213, 'n_r': 1.42, 'n_i': 0.007, 'n_scale': 0.0015},
                          {'model': 'urban', 'mode': 'coarse', 'rv': 3.3252, 'rv_scale': 0.1411, 'sig': 0.7595, 'sig_scale': 0.1638, 'Cv': 0.0934, 'Cv_scale': 0.6394, 'n_r': 1.42, 'n_i': 0.007, 'n_scale': 0.0015},
                          {'model': 'smoke', 'mode': 'fine', 'rv': 0.1335, 'rv_scale': 0.0096, 'sig': 0.3834, 'sig_scale': 0.0794, 'Cv': 0.1748, 'Cv_scale': 0.8914, 'n_r': 1.51,'n_i': 0.02,   'n_scale': 0},
                          {'model': 'smoke', 'mode': 'coarse', 'rv': 3.4479, 'rv_scale': 0.9489, 'sig': 0.7433, 'sig_scale': 0.0409, 'Cv': 0.1043, 'Cv_scale': 0.6824, 'n_r': 1.51, 'n_i': 0.02, 'n_scale': 0},
                          {'model': 'dust', 'mode': 'fine', 'rv': 0.1416, 'rv_scale': -0.0519, 'sig': 0.7561, 'sig_scale': 0.148, 'Cv': 0.087, 'Cv_scale': 1.026, 'n_r': 1.48, 'n_i': 0.0025, 'n_scale': (-0.021,0.132)},
                          {'model': 'dust', 'mode': 'coarse', 'rv':2.2, 'rv_scale': 0, 'sig': 0.554, 'sig_scale': -0.0519, 'Cv': 0.6786, 'Cv_scale': 1.0569, 'n_r': 1.48, 'n_i': 0.0025, 'n_scale': (-0.021,0.132)},
                         ])
        self.model_parameters = models
        self.aod = aod   
        
        r_range = np.array(diameter_range[:2])/2/1e3
        r = np.logspace(np.log10(r_range[0]), np.log10(r_range[1]), diameter_range[2])
                
        
        #####
        bins, names = db.bincenters2binsANDnames(r*2*1e3)
        dists = {}
        for mo in models.model.unique():
            mos = models[models.model == mo]
        
            dist = np.zeros(r.shape)
            for idx,row in mos.iterrows():
                if row.model == 'dust':
                    rv = row.rv * aod**row.rv_scale
                    sig = row.sig  * aod**row.sig_scale
                else:
                    rv = row.rv + (row.rv_scale * aod)
                    sig = row.sig + (row.sig_scale * aod)
        
                Cv = row.Cv * aod**row.Cv_scale
        
                dist += Cv/(np.sqrt(2* np.pi) * sig) * np.exp(- (np.log(r) - np.log(rv))**2/(2 * sig**2))
            dist = sd.SizeDist(pd.DataFrame([dist], columns=names), bins, 'dVdlogDp')
            dists[mo] = dist
            
        self.models = dists
        for mo in dists:
            dist = dists[mo]
            setattr(self, mo, dists[mo])
            
            
    def plot_all(self, ax = None, moment = 'V'):
        """
        plot models

        Parameters
        ----------
        ax : TYPE, optional
            For plotting on existing axis. The default is None.
        moment : TYPE, optional
            For plotting in (V)olume, (S)urface, or (N)umber moment; all in 
            dx/dlogdp . The default is 'V'.

        Returns
        -------
        a : TYPE
            DESCRIPTION.

        """
        if isinstance(ax, type(None)):
            f,a = plt.subplots()
        else:
            a = ax
        
        dists = self.models
        for mo in dists:
            dist = dists[mo]
            if moment == 'N':
                dist = dist.convert2dNdlogDp()
            elif moment == 'S':
                dist = dist.convert2dSdlogDp()
            elif moment == 'V':
                dist = dist.convert2dVdlogDp()
            else:
                raise ValueError(f'Kwarg "moment" has to be "N", "S", or "V", not "{moment}"')
            # g, = a.plot(r*2, dist)
            f,a = dist.plot(ax = a)
            g = a.get_lines()[-1]
            g.set_label(mo)
            
        # a.set_xscale('log')
        a.legend()
        return a
        