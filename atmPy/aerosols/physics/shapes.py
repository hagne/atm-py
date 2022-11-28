# import pkg_resources
import xarray as _xr
# import importlib.resources
# import pkgutil
# from io import StringIO
# import os
import pathlib as pl
import numpy as np
import xarray as xr
import scipy as sp


def _load_davies1979():
    """
    This dataset is an extraction of data from 
    
    Davies, C. N. (1979). Particle-fluid interaction. Journal of Aerosol 
    Science, 10(5), 477–513. https://doi.org/10.1016/0021-8502(79)90006-5
    
    For details on how the dataset extracted/generated check the notebook
    fundgrube/aerosols/sphericity/aspectratio_vs_dynamicshapefactor.ipynb

    Returns
    -------
    out : TYPE
        DESCRIPTION.

    """
    p2f = pl.Path(__file__).parent.joinpath('data/spheroid2dynamic_shape_Davies1979.nc')
    out = _xr.open_dataset(p2f)
    return out

class AspectRatio2DynamicShapeFactor():
    def __init__(self,basedon="davies1979"):
        self.basedon = basedon
        self._dataset = None
        return
    
    @property
    def dataset(self):
        if isinstance(self._dataset, type(None)):
            assert(self.basedon == 'davies1979')
            self._dataset = _load_davies1979()
        return self._dataset
    
    def shape_distribution2effective_dsf(self, limits, distribution_function = 'squared', 
                                         # scale = 5,
                                         ):
        """
        Determines the effective dynamic shape factor from a distribution of
        aspect ratios based on the shape parameter.

        Parameters
        ----------
        limits : TYPE
            DESCRIPTION.
        distribution_function : TYPE, optional
            DESCRIPTION. The default is 'squared'.

        Returns
        -------
        dsf_eff : TYPE
            DESCRIPTION.

        """
        # testing if sampling rate is uniform
        s2d = self.dataset
        nda = s2d.shape_parameter.values
        spacing = nda[1:] - nda[:-1]
        assert(np.std(spacing)/np.mean(spacing) < 1e-3), 'Data is not equally distributed, this will case a bias!! You probably changed the data at some point, because that did not happen with the original data!!!'
        # limits
        # if limits array is only one element long (happens in optimization routine)
        if hasattr(limits, '__iter__'):
            if len(limits) < 2:
                limits = limits[0]
         
        # if limits is single number make limits symetric around 0
        if not hasattr(limits, '__iter__'):
            limits = [-limits,limits]
        
        
        assert(limits[0] > float(self.dataset.shape_parameter.min())), 'limits[0] out of bounds'
        assert(limits[1] < float(self.dataset.shape_parameter.max())), 'limits[1] out of bounds'
        self.tp_limits = limits
        
            
        
        # effective dynamicshapefactor
        if distribution_function == 'squared':
            # nofv = (limits[1] - limits[0])/np.mean(spacing) * scale
            #int(np.ceil(nofv)) 
            # self.tp_nsp = int(np.ceil(nofv))
            
            noofv = 250 # the number of 250 results in an error smaller than 0.5% for dsf between 1.05, and 1.5
            new_shape_params = np.linspace(limits[0], limits[1], noofv)
            new_dsf = self.dataset.dynamic_shape_factor.interp(shape_parameter=new_shape_params, method = 'quadratic')
            self.tp_sp_dist = new_dsf
            dsf_eff = float(new_dsf.sum()/new_dsf.shape[0])
            
        elif distribution_function == 'squared1':
            sp_dist = np.ones(s2d.shape_parameter.shape)
            sp_dist[np.where(s2d.shape_parameter < limits[0])] = 0
            sp_dist[np.where(s2d.shape_parameter > limits[1])] = 0
            
            sp_dist = xr.DataArray(sp_dist, 
                    coords=[s2d.shape_parameter,],
                    )
            
            self.tp_sp_dist = sp_dist
            dsf_eff = float((s2d.dynamic_shape_factor * sp_dist).sum()/ sp_dist.sum())
            
        else:
            assert(False), 'not implemented yet'
            
        return dsf_eff
            
    def effective_dsf2shape_distribution(self, dsf, bounds=(1, 10)):
        res = sp.optimize.least_squares(lambda sp: (self.shape_distribution2effective_dsf(sp) - dsf) * 1e2, [5,], bounds=bounds) 
        return float(res.x)

# Sphereoids
def volume_equivalent_radius_spheroid(equatorial_axis, polar_axis):
    a = equatorial_axis
    c = polar_axis
    return  (a**2 * c)**(1/3)

def projection_area_equivalent_radius_spheroid(axis1, axis2):
    a = axis1
    c = axis2
    return (a * c) ** (1 / 2)


