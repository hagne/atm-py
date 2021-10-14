#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 21 10:20:32 2021

@author: hagen
"""

import pandas as pd
import pathlib as pl
import numpy as np
    
def fixtime(df):
    if 'EPOCH' in df.columns:
        dt = pd.to_datetime(df.EPOCH, unit = 's')
        col = 'EPOCH'
    elif 'DateTimeUTC' in df.columns:
        dt = pd.to_datetime(df.DateTimeUTC)
        col = 'DateTimeUTC'
    df.index = dt
    df.index.name = 'DateTimeUTC'
    df.drop(col, axis = 1, inplace = True)

def generate_cpd3_export_command(num_days = None, 
                                 start = None, 
                                 end = None, #'2021-05-13T00:00:00'
                                 status = 'clean',
                                 p2f_pops = None,#'/home/grad/htelg/cpd3_exports/cpd3_qc_qa_pops.csv',
                                 p2f_smps = None,#'/home/grad/htelg/cpd3_exports/cpd3_qc_qa_smps.csv',
                                 p2f_neph = None,#'/home/grad/htelg/cpd3_exports/cpd3_qc_qa_neph.csv',
                                ):
    """
    Generate the command that needs to be executed in the commandline on aero

    Parameters
    ----------
    num_days : TYPE, optional
        DESCRIPTION. The default is None.
    start : TYPE, optional
        DESCRIPTION. The default is None.
    end : TYPE, optional
        DESCRIPTION. The default is None.
    #'2021-05-13T00 : 00:00'                                 
    status : 'str' ['clean', 'raw'], optional
        Currently this effects only the POPS data (might not be true anymore). If raw, the data is not 
        corrected (e.g. for T,P) and qa/qc has not been applied The default is 'clean'.
    p2f_pops : TYPE, optional
        DESCRIPTION. The default is None.
    #'/home/grad/htelg/cpd3_exports/cpd3_qc_qa_pops.csv' : TYPE
        DESCRIPTION.
    p2f_smps : TYPE, optional
        DESCRIPTION. The default is None.
    #'/home/grad/htelg/cpd3_exports/cpd3_qc_qa_smps.csv' : TYPE
        DESCRIPTION.
    p2f_neph : TYPE, optional
        DESCRIPTION. The default is None.
    #'/home/grad/htelg/cpd3_exports/cpd3_qc_qa_neph.csv' : TYPE
        DESCRIPTION.
     : TYPE
        DESCRIPTION.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    end_dt = None
    start_dt = None
    num_days_td = None
    if not isinstance(end, type(None)):
        end_dt = pd.to_datetime(end)
    if not isinstance(start, type(None)):
        start_dt = pd.to_datetime(start)
    if not isinstance(num_days, type(None)):
        num_days_td = pd.to_timedelta(num_days,'d')
        
    if np.all([not isinstance(end, type(None)), not isinstance(start, type(None)), not isinstance(num_days, type(None))]):
        assert(False), "one of [num_days, start, end] needs to be None"
    elif np.all([isinstance(end, type(None)), isinstance(start, type(None)), isinstance(num_days, type(None))]):
        assert(False), "at least on of of [num_days, start, end] needs to be given."
    elif np.all([not isinstance(end, type(None)), not isinstance(start, type(None))]):
        pass
    elif np.all([not isinstance(end, type(None)), not isinstance(num_days, type(None))]):
        start_dt = end_dt - num_days_td
    elif np.all([not isinstance(start, type(None)), not isinstance(num_days, type(None))]):
        end_dt = start_dt + num_days_td
    elif not isinstance(num_days, type(None)):
        assert(np.all([isinstance(end, type(None)), isinstance(start, type(None))])), 'end and start are not None, this should not happen'
        end_dt = pd.Timestamp.now()
        start_dt = end_dt - num_days_td
    elif not isinstance(start, type(None)):
        end_dt = pd.Timestamp.now()
    else:
        assert(False), "not possible"
        
    # if isinstance(end, type(None)):
    #     assert(np.all([isinstance(end, type(None)), isinstance(start, type(None)))), 'end and start are not None, this should not happen'
    #     end_dt = pd.Timestamp.now()
    #     start_dt = end_dt - num_days_td
    # else:
    #     now = pd.to_datetime(now)
        
    # yd= now - pd.to_timedelta(num_days,'d')
    start = f'{start_dt.year}-{start_dt.month:02d}-{start_dt.day:02d}T00:00:00Z'
    end = f'{end_dt.year}-{end_dt.month:02d}-{end_dt.day:02d}T00:00:00Z'  
        
    coms = []
    if not isinstance(p2f_pops, type(None)):
        p2f_pops = pl.Path(p2f_pops)
        coms.append(f"da.export bos '.*_N21' {start} {end} {status}  > {p2f_pops}")
    if not isinstance(p2f_smps, type(None)):
        p2f_smps = pl.Path(p2f_smps)
        coms.append(f"da.export bos '.*_N11' {start} {end} {status} > {p2f_smps}")
    if not isinstance(p2f_smps, type(None)):
        p2f_neph = pl.Path(p2f_neph)
        coms.append(f"da.get bos S11a {start} {end} {status} | da.avg --interval=5m | da.export --mode=R --fill=1h > {p2f_neph}")


    # com_pops = f"da.export bos '.*_N21' {start} {end}  > {p2f_pops}"
    # com_smps = f"da.export bos '.*_N11' {start} {end}  > {p2f_smps}"
    # com_neph = f"da.get bos S11a {start} {end} raw | da.avg --interval=5m | da.export --mode=R --fill=1h > {p2f_neph}"
    # coms = [com_pops, com_smps, com_neph]
    return '; '.join(coms)