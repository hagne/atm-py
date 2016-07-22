import numpy as np
import pandas as pd
test_data_folder = './test_data/'
#### data archives
######## ARM
from atmPy.data_archives.arm import read_data

def test_1twr10xC1():
    out = read_data.read_cdf(test_data_folder, data_product='1twr10xC1')
    out = out['1twr10xC1']

    # rh
    soll = pd.read_csv(test_data_folder + '1twr10xC1_rh.csv', index_col=0,
                       dtype={'rh_25m': np.float32, 'rh_60m': np.float32}
                       )
    assert np.all(out.relative_humidity.data == soll)

    # temp
    soll = pd.read_csv(test_data_folder + '1twr10xC1_temp.csv', index_col=0,
                       dtype={'temp_25m': np.float32, 'temp_60m': np.float32}
                       )
    assert np.all(out.temperature.data == soll)

    # vapor pressure
    soll = pd.read_csv(test_data_folder + '1twr10xC1_p_vapor.csv', index_col=0,
                       dtype={'vap_pres_25m': np.float32, 'vap_pres_60m': np.float32}
                       )

    assert np.all(out.vapor_pressure.data == soll)