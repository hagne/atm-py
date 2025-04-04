# This module was written under the assumption that the arm data is extremely uniform ... whell its not, so its bett
# to write an opener for each separately
#
# from atmPy.data_archives.arm._netCDF import ArmDataset as _Dataset
import os as _os
# from atmPy.data_archives.arm.file_io.products import _tdmasize,_tdmaapssize,_tdmahyg,_aosacsm, _noaaaos, _1twr10xC1, _aipfitrh1ogrenC1, aosaps, aossmps
import pandas as _pd
import pylab as _plt
import warnings
# from atmPy.data_archives.arm import tools
import pathlib as _pl
import numpy as _np
import xarray as _xr
# import pdb as _pdb

arm_products = {}
# arm_products = {'aosaps':     {'module': aosaps},
#                 'aossmps':     {'module': aossmps},
#                 'tdmasize':   {'module': _tdmasize},
#                 'tdmaapssize':{'module': _tdmaapssize},
#                 'tdmahyg':    {'module': _tdmahyg},
#                 'aosacsm':    {'module': _aosacsm},
#                 'noaaaos':    {'module': _noaaaos},
#                 '1twr10xC1':  {'module': _1twr10xC1},
#                 'aipfitrh1ogrenC1': {'module': _aipfitrh1ogrenC1}
#                 }

class Generic_Arm_Product(object):
    def __init__(self, dataset):
        self.dataset = dataset
        pass
    
class Generic_Arm_Variable(object):
    def __init__(self, parent, data, qc, instance_of_variable = None):
        self.product = parent
        self._orig_data = data
        self.qc = qc
        self._data_good = None
        self._data_intermediate = None
        self._data_intermediate_good= None
        self._data_intermediate= None
        self._data_bad= None
        self._data_all = None
        
        if isinstance(instance_of_variable, type(None)):
            self.instancifyer = lambda x: x
        else:
            self.instancifyer = lambda x: instance_of_variable(x)

    @property
    def data_all(self):
        if isinstance(self._data_all, type(None)):
            self._data_all = self.apply_mask([])
        return self._data_all
    
    @property
    def data_good(self):
        if isinstance(self._data_good, type(None)):
            self._data_good = self.apply_mask(['Bad', 'Indeterminate'])
        return self._data_good
    
    @property
    def data_intermediate_good(self):
        if isinstance(self._data_intermediate_good, type(None)):
            self._data_intermediate_good = self.apply_mask(['Bad',])
        return self._data_intermediate_good
    
    @property
    def data_intermediate(self):
        if isinstance(self._data_intermediate, type(None)):
            self._data_intermediate = self.apply_mask(['Bad', 'Good'])
        return self._data_intermediate
    
    @property
    def data_bad(self):
        if isinstance(self._data_bad, type(None)):
            self._data_bad = self.apply_mask(['Good', 'Indeterminate'])
        return self._data_bad
    
    def apply_mask(self, dataqualit2mask):
        data_masked = self._orig_data.copy()
        for bd in dataqualit2mask:
            data_masked[self.qc == bd] = _np.nan
        return self.instancifyer(data_masked)
    
class Mfrsr7nchaod1michC1_AOD(Generic_Arm_Variable):
    def __init__(self,*args, **kwargs):
        import atmPy.aerosols.physics.column_optical_properties as atmcop
        kwargs['instance_of_variable'] = atmcop.AOD_AOT
        super().__init__(*args, **kwargs)
        ds = self.product.dataset
        wavelength_nominal_nm = {k.split('_')[0]: int(ds.attrs[k].split()[0]) for k in ds.attrs if 'CWL_nominal' in k}
        wavelength_measured_nm = {k.split('_')[0]: float(ds.attrs[k].split()[0]) for k in ds.attrs if 'CWL_measured' in k}

        self._orig_data.rename(wavelength_nominal_nm, axis = 1, inplace = True)
        self.qc.rename(wavelength_nominal_nm, axis = 1, inplace = True)

def generic_arm_file_reader(path2file, variables = [], variable_names = [], 
                            variable_instance = Generic_Arm_Variable,
                            mask_data = ['Bad', 'Indeterminate'],
                            verbose = True):

    p2f = _pl.Path(path2file)

    ds = _xr.open_dataset(p2f)
    
    if ds.attrs['platform_id'] == 'mfrsr7nchaod1mich':
        if verbose:
            print('product is mfrsr7nchaod1mich')
        variable_instance = Mfrsr7nchaod1michC1_AOD
        variables = ['aerosol_optical_depth_filter']
        variable_names = ['AOD']
    else:
        variable_instance = Generic_Arm_Variable
        variables = None
        
        
    prod = Generic_Arm_Product(ds)
    
    if isinstance(variables, str):
        variables = [variables,]
    
    for e,variable in enumerate(variables):
        aod_vars = [v for v in list(ds.variables) if variable in v]
        aod_vars_qc = [v for v in aod_vars if 'qc_' in v]
        aod_vars = [v for v in aod_vars if 'qc_' not in v]

        aod_df = ds[aod_vars].to_pandas()
        aod_df_qc = ds[aod_vars_qc].to_pandas()

        aod_df = _pd.DataFrame()
        aod_assess_df = _pd.DataFrame()

        for var in aod_vars:
            # break

            # column name simplification
            col = var.split('_')[-1]

            # variable
            aod_df[col] = ds[var].to_pandas()

            # qc
            qc = ds[f'qc_{var}'].to_pandas()

            # how many bits are there?
            noofbits = _np.max([int(k.split('_')[1]) for k in ds[f'qc_{var}'].attrs if 'bit' in k])

            #qc to bits
            qcb = qc.apply(lambda q: _np.binary_repr(q,width=noofbits))

            # get the bit assessment (good, intermediate, bad)
            bitassemsments = [ds[f'qc_{var}'].attrs[k] for k in ds[f'qc_{var}'].attrs if ('bit' in k) and ('assessment' in k)]
            qca = qcb.apply(lambda val: sorted([bitassemsments[e] if qa == '1' else 'Good' for e,qa in enumerate(val[::-1])])[0])
            aod_assess_df[col] = qca



        var  = variable_instance(prod, aod_df, aod_assess_df)
        setattr(prod, variable_names[e], var)
        
    return prod



#### older stuff
def open_path(path, **kwargs):
    """

    Parameters
    ----------
    path: str
    window: tuple
        e.g. ('1990-01-01','2030-01-01')
    kwargs: all other kwargs the particular file might take, see the module for details

    Returns
    -------

    """
    info = _tools.path2info(path)
    module = arm_products[info['product']]['module']
    out = module.open_path(path, **kwargs)
    return out

def check_availability(folder,
                       data_product = None,
                       site = 'sgp',
                       time_window = ('1990-01-01','2030-01-01'),
                       custom_product_keys = False,
                       ignore_unknown = True,
                       verbose = False):

    fname = _os.listdir(folder)
    index = _pd.date_range('1990-01-01','2030-01-01', freq = 'D')
    df = _pd.DataFrame(index = index)

    for f in fname:
        if verbose:
            print('\n', f)

        # error handling: test for netCDF file format
        if _os.path.splitext(f)[-1] != '.cdf':
            txt = '\t %s is not a netCDF file ... skipping'%f
            if verbose:
                print(txt)
            continue

        site_check = _is_site(f,site,verbose)
        if not site_check:
            continue

        date = _is_in_time_window(f,time_window,verbose)
        if not date:
            continue

        product_id = _is_in_product_keys(f, ignore_unknown, verbose, custom_product_keys = custom_product_keys)
        if not product_id:
            continue

        if not _is_desired_product(product_id,data_product,verbose):
            continue

        if product_id not in df.columns:
            df[product_id] = _pd.Series(1, index = [date])
        else:
            df[product_id][date] = 1

    df = df.sort_index(axis=1)

    for e,col in enumerate(df.columns):
        df[col].values[df[col].values == 1] = e+1


    f,a = _plt.subplots()
    for col in df.columns:
        a.plot(df.index,df[col], lw = 35, color = [0,0,1,0.3])

    a.set_ylim((0.1,df.shape[1] + 0.9))
    bla = range(1,df.shape[1]+1)
    a.yaxis.set_ticks(bla)
    a.yaxis.set_ticklabels(df.columns)

    f.autofmt_xdate()

    f.tight_layout()
    return df, a


def read_cdf(fname,
             site = 'sgp',
             data_product = None,
             time_window = None,
             data_quality = 'good',
             data_quality_flag_max = None,
             check_in_subfolder = False,
             concat = True,
             ignore_unknown = False,
             leave_cdf_open = False,
             verbose = False,
             error_bad_file = True,
             ):
    """
    Reads ARM NetCDF file(s) and returns a containers with the results.

    Parameters
    ----------
    fname: str or list of str.
        Either a file, directory, or list of files. If directory name is given
        all files in the directory will be considered.
    data_product: str.
        To see a list of allowed products look at the variable arm_products.
    time_window: tuple of str.
        e.g. ('2016-01-25 15:22:40','2016-01-29 15:00:00').
        Currently the entire day is considered, no need to use exact times.
    check_in_subfolder: bool
        if to check for files in subfolder. If True only files in subfolders are considered and only one level deep!
    concat
    ignore_unknown
    verbose

    Returns
    -------

    """

    # list or single file
    # if type(fname) == str:
    #     if _os.path.isdir(fname):
    #         f = _os.listdir(fname)
    #         fname = [fname + i for i in f]
    #     else:
    #         fname = [fname]

    if type(fname) == str:
        if _os.path.isdir(fname):

            if check_in_subfolder:
                fnames = []
                for fol in _os.listdir(fname):
                    if _os.path.isdir(fname + fol):
                        f = _os.listdir(fname + fol)
                        fnames += [fname + fol + '/' + i for i in f]
                fname = fnames
            else:
                f = _os.listdir(fname)
                fname = [fname + i for i in f]
        else:
            fname = [fname]

    if len(fname) > 1 and leave_cdf_open:
        txt = "leave_cdf_open can only be true if the number of files is one ... leave_cdf_open = False"
        warnings.warn(txt)
        leave_cdf_open = False

    if type(data_product) == str:
        data_product = [data_product]
    products = {}

    #loop thru files
    no_valid = 0
    for f in fname:
        if verbose:
            print('\n', f)

        # error handling: test for netCDF file format
        # _pdb.set_trace()
        if _os.path.splitext(f)[-1] != '.cdf':
            txt = '\t %s is not a netCDF file ... skipping'%f
            if verbose:
                print(txt)
            continue

        if not _is_in_time_window(f,time_window,verbose):
            continue

        product_id = _is_in_product_keys(f, ignore_unknown, verbose)
        if not product_id:
            continue

        site_check = _is_site(f,site,verbose)
        if not site_check:
            continue

        if not _is_desired_product(product_id,data_product,verbose):
            continue

        if product_id not in products.keys():
            products[product_id] = []


        arm_file_object = arm_products[product_id]['module'].ArmDatasetSub(f,
                                                                           data_quality = data_quality,
                                                                           data_quality_flag_max = data_quality_flag_max,
                                                                           error_bad_file = error_bad_file)

        if not leave_cdf_open:
            arm_file_object._close()

        # if there was an error in reading the time stamp, the file will be discarded
        # import pdb
        # pdb.set_trace()

        no_valid += 1
        if arm_file_object._parsing_error:
            continue

        products[product_id].append(arm_file_object)

    if len(fname) == 1:
        if not no_valid:
            txt = '%s is either not the right file format or does not fall into the enquiry specifications'%(fname[0])
            raise ValueError(txt)
        return arm_file_object

    else:
        if concat:
            for pf in products.keys():
                products[pf] = arm_products[pf]['module']._concat_rules(products[pf])
        return products


def _is_desired_product(product_id, data_product, verbose):
    out = True
    if data_product:
        if product_id not in data_product:
            if verbose:
                print('Not the desired data product ... skip')
            out = False
    return out

def _is_site(f,site,verbose):
    out = True
    fnt = _os.path.split(f)[-1].split('.')[0]
    site_is = fnt[:3]
    if site:
        if site_is != site:
            out = False
            if verbose:
                txt = 'Has wrong site_id (%s) ... skip!'%(site_is)
                print(txt)
    return out

def _is_in_product_keys(f, ignore_unknown,verbose, custom_product_keys = False):

    fnt = _os.path.split(f)[-1].split('.')
    product_id = False
    for prod in arm_products.keys():
        if prod in fnt[0]:
            product_id = prod
            break

    if custom_product_keys:
        for prod in custom_product_keys:
            if prod in fnt[0]:
                product_id = prod
                return product_id

    if not product_id:
        txt = '\t has no ncattr named platform_id. Guess from file name failed ... skip'
        if verbose:
            print(txt)
    else:
        if product_id not in arm_products.keys():
            txt = 'Platform id %s is unknown.'%product_id
            product_id = False
            if ignore_unknown:
                if verbose:
                    print(txt + '... skipping')
            else:
                raise KeyError(txt)
    return product_id

def _is_in_time_window(f,time_window, verbose):
    out = True
    if time_window:
        fnt = _os.path.split(f)[-1].split('.')
        ts = fnt[-3]
        file_start_data = _pd.to_datetime(ts)
        start_time = _pd.to_datetime(time_window[0])
        end_time = _pd.to_datetime(time_window[1])
        dt_start = file_start_data - start_time
        dt_end = file_start_data - end_time
        out = file_start_data

        if dt_start.total_seconds() < -86399:
            if verbose:
                print('outside (before) the time window ... skip')
            out = False
        elif dt_end.total_seconds() > 86399:
            if verbose:
                print('outside (after) the time window ... skip')
            out = False
    return out