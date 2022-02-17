from io import StringIO as _StringIO
import pandas as _pd
from atmPy.aerosols.instruments.nephelometer import data
import matplotlib.pylab as _plt

class Nephelometer(object):
    def __init__(self):
        #         self._absorbtion_coeff = None
        self._scattering_coeff = None
        self._hemisphericbackscatt_coeff = None
        self._RH_in_instrument = None
        self.channels = Channels()



    @property
    def scattering_coeff(self):
        return self._scattering_coeff

    @scattering_coeff.setter
    def scattering_coeff(self, value):
        self._scattering_coeff = value

    @property
    def hemisphericbackscatt_coeff(self):
        return self._hemisphericbackscatt_coeff

    @hemisphericbackscatt_coeff.setter
    def hemisphericbackscatt_coeff(self, value):
        self._hemisphericbackscatt_coeff = value

    @property
    def RH_in_instrument(self):
        return self._RH_in_instrument

    @RH_in_instrument.setter
    def RH_in_instrument(self, value):
        self._RH_in_instrument = value

class TandemNephelometer(object):
    def __init__(self):
        self._nephelometer_dry = None
        self._nephelometer_wet = None

class Channels(object):
    def __init__(self):
        buff = _StringIO(data.data_450)
        df_450 = _pd.read_csv(buff, sep='  ', index_col=1, names=['450', 'wavelength'], engine='python')
        df_450.sort_index(inplace=True)
        self.channel_450 = df_450

        buff = _StringIO(data.data_550)
        df_550 = _pd.read_csv(buff, sep='  ', index_col=1, names=['550', 'wavelength'], engine='python')
        df_550.sort_index(inplace=True)
        self.channel_550 = df_550

        buff = _StringIO(data.data_700)
        df_700 = _pd.read_csv(buff, sep='  ', index_col=1, names=['700', 'wavelength'], engine='python')
        df_700.sort_index(inplace=True)
        self.channel_700 = df_700

    def plot(self, ax = None):
        if not ax:
            f,a = _plt.subplots()
        self.channel_450.plot(ax = a)
        self.channel_550.plot(ax = a)
        self.channel_700.plot(ax = a)


