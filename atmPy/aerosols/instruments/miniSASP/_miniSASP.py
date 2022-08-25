from atmPy.general import timeseries
from atmPy.radiation import solar
import pandas as pd
from atmPy.radiation.rayleigh import bucholtz_rayleigh as _bray
# from scipy import signal
import numpy as np
# import matplotlib.pylab as plt
# from atmPy.tools import plt_tools

class MiniSASP(timeseries.TimeSeries):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._data_period = 0.03
        self.wavelength_channels = {'PhotoA': '550.4',
                                    'PhotoB': '460.3',
                                    'PhotoC': '671.2',
                                    'PhotoD': '860.7'}

        self.__sun_intensities = None

        self.__sun_elevation = None

        self.__housekeeping = None

        self.__od_amf = None

        self.__air_mass_factor = None

        self.__od = None

        self.optical_depth_rayleigh_offsets =  {'550.4': 0,
                                                '460.3': 0,
                                                '671.2': 0,
                                                '860.7': 0}

        self.__od_ray_orig = None
        self.__od_ray_offset_last = None
        self.optical_depth_amf_offsets = {'460.3': 0, '550.4': 0, '671.2': 0, '860.7': 0}
        self.__od_afm_offset_last = None
        self.__aod = None
        self.__od_amf_orig = None



    # @property
    # def sun_intensities(self, which='short', min_snr=10, moving_max_window=23):
    #     """ Finds the peaks in all four photo channels (short exposure). It also returns a moving maximum as guide to
    #      the eye for the more "real" values.
    #
    #      Parameters
    #      ----------
    #      which: 'long' or 'short'.
    #         If e.g. PhotoA (long) or PhotoAsh(short) are used.
    #      min_snr: int, optional.
    #         Minimum signal to noise ratio.
    #      moving_max_window: in, optionl.
    #         Window width for the moving maximum.
    #
    #
    #      Returns
    #      -------
    #      TimeSeries instance (AtmPy)
    #     """
    #     if not self.__sun_intensities:
    #         def _moving_max(ds, window=3):
    #             out = pd.DataFrame(ds, index=ds.index)
    #             out = pd.rolling_max(out, window)
    #             out = pd.rolling_mean(out, int(window / 5), center=True)
    #             return out
    #         moving_max_window = int(moving_max_window / 2.)
    #         # till = 10000
    #         # photos = [self.data.PhotoAsh[:till], self.data.PhotoBsh[:till], self.data.PhotoCsh[:till], self.data.PhotoDsh[:till]]
    #         if which == 'short':
    #             photos = [self.data.PhotoAsh, self.data.PhotoBsh, self.data.PhotoCsh, self.data.PhotoDsh]
    #         elif which == 'long':
    #             photos = [self.data.PhotoA, self.data.PhotoB, self.data.PhotoC, self.data.PhotoD]
    #         else:
    #             raise ValueError('which must be "long" or "short" and not %s' % which)
    #         # channels = [ 550.4, 460.3, 860.7, 671.2]
    #         df_list = []
    #         for e, data in enumerate(photos):
    #             pos_indx = signal.find_peaks_cwt(data, np.array([10]), min_snr=min_snr)
    #             out = pd.DataFrame()
    #             # return out
    #             # break
    #             out[str(self.wavelength_channels[e])] = pd.Series(data.values[pos_indx], index=data.index.values[pos_indx])
    #
    #             out['%s max' % self.wavelength_channels[e]] = _moving_max(out[str(self.wavelength_channels[e])], window=moving_max_window)
    #             df_list.append(out)
    #
    #         out = pd.concat(df_list).sort_index()
    #         # out = out.groupby(out.index).mean()
    #         self.__sun_intensities = Sun_Intensities_TS(out)
    #     return self.__sun_intensities

    @property
    def sun_intensities(self):
        threshold_scale = 5
        verbose = False
        if not self.__sun_intensities:
            start_time = self.get_timespan()[0].to_datetime64()
            # channel_list = ['PhotoAsh', 'PhotoBsh', 'PhotoCsh', 'PhotoDsh']
            channel_list = [i + 'sh' for i in np.array(sorted(self.wavelength_channels.items(), key=lambda x: x[1]))[:, 0]]
            wl_list = np.array(sorted(self.wavelength_channels.items(), key=lambda x: x[1]))[:, 1]
            sun_max_list = []
            for channel in channel_list:
                dft = self.drop(channel, inverse = True)
                #         dft = dft.zoom_time(start = '2015-04-21 08:01:00.0000',end = '2015-04-21 08:04:00.0000')
                dft = dft.data

                # threshold = 50 * np.median(dft)
                threshold = threshold_scale * np.median(dft.dropna())
                if verbose:
                    print('threshold: %s' % threshold)

                above_threshold = dft > threshold
                above_threshold_inv = above_threshold != True

                # mean operation can not be performed on time, therefore I convert the time into into milliseconds since some arbitrary date
                dft['time_since_start'] = (dft.index - start_time) / np.timedelta64(1, 'ms')

                # group data of sections above shreshold (... the peak)
                atic = above_threshold_inv.cumsum()
                atic[above_threshold_inv == True] = 0
                dft['above_threshold_inv'] = atic
                grp = dft.groupby('above_threshold_inv')

                # mean of groups provides center time
                grp_mean = grp.mean()
                grp_mean = grp_mean.drop([channel], axis=1)

                # max gives the max value in each group
                grp_max = grp.max()
                grp_max = grp_max.drop(['time_since_start'], axis=1)

                # Merge the two in to one time series ... including retriveing the time from the "time since ..."
                sun_ints = pd.DataFrame(index=grp_mean.time_since_start)
                sun_ints['max_ints'] = grp_max.values
                # # sun_ints['time_since_start'] = grp_mean
                sun_ints.index = (sun_ints.index.values * np.timedelta64(1, 'ms')) + start_time
                sun_ints.drop([sun_ints.index[0]], inplace=True)
                sun_max_list.append(sun_ints)

            sun_ints_ts = sun_max_list[0]
            sun_ints_ts.columns = [wl_list[0]]#[channel_list[0]]

            for e, sm in enumerate(sun_max_list[1:]):
                sun_ints_ts[wl_list[e + 1]] = sm.reindex(sun_ints_ts.index, method='nearest')

            self.__sun_intensities = Sun_Intensities_TS(sun_ints_ts)
        return self.__sun_intensities


    @sun_intensities.setter
    def sun_intensities(self,data):
        self.__sun_intensities = data


    @property
    def optical_depth(self):
        if not self.__optical_depth:
            pass


    @property
    def housekeeping(self):
        """Although the position of the instrument is recorded by the internal instruments own GPS, precision is low and altitude is missing. Provide the information by setting this property.
        alignment will be taken care of."""
        return self.__housekeeping

    @housekeeping.setter
    def housekeeping(self, ts):
        # gps_cols = ['Lat', 'Lon', 'Altitude']
        # ts = ts._del_all_columns_but(gps_cols)
        # cols = ts.data.columns
        # merged = self.merge(ts, recognize_gaps=False)
        # self.__housekeeping = merged._del_all_columns_but(cols)
        self.__housekeeping = ts



    @property
    def sun_elevetion(self):
        """
        doc is not correct!!!

        This function uses telemetry data from the airplain (any timeseries including Lat and Lon) to calculate
        the sun's elevation. Based on the sun's elevation an airmass factor is calculated which the data is corrected for.

        Arguments
        ---------
        sun_intensities: Sun_Intensities_TS instance
        picco: any timeseries instance containing Lat and Lon
        """
        if not self.housekeeping:
            txt = 'For this calculation we need information on Lat, Lon, Altitude. Please set the attribute housekeeping with a timeseries that has these informations'
            raise AttributeError(txt)

        if not self.__sun_elevation:
            cols = self.housekeeping.data.columns
            merged = self.merge(self.housekeeping, recognize_gaps=False)
            merged = merged._del_all_columns_but(cols)

            ts = solar.get_sun_position_TS(merged)
            ts = ts._del_all_columns_but(['Solar_position_elevation'])
            self.__sun_elevation = ts
        # picco_t = timeseries.TimeSeries(picco.data.loc[:, ['Lat', 'Lon', 'Altitude']])  # only Altitude, Lat and Lon
        # sun_int_su = self.merge(picco_t)
        # out = sun_int_su.get_sun_position()
        # #     sun_int_su = sun_int_su.zoom_time(spiral_up_start, spiral_up_end)
        # arrays = np.array([sun_int_su.data.index, sun_int_su.data.Altitude, sun_int_su.data.Solar_position_elevation])
        # tuples = list(zip(*arrays))
        # index = pd.MultiIndex.from_tuples(tuples, names=['Time', 'Altitude', 'Sunelevation'])
        # sun_int_su.data.index = index
        # sun_int_su.data = sun_int_su.data.drop(
        #     ['Altitude', 'Solar_position_elevation', 'Solar_position_azimuth', 'Lon', 'Lat'], axis=1)
        return self.__sun_elevation


    @property
    def air_mass_factor(self):
        if not self.__air_mass_factor:
            self.__air_mass_factor = timeseries.TimeSeries(1. / np.sin(self.sun_elevetion.data))
            self.__air_mass_factor._data_period = self.sun_elevetion._data_period
        return self.__air_mass_factor


    @property
    def optical_depth_amf(self):
        """OD * airmassfactor + unkonwn offset. after determining the offset you might want to set the property sup_offsets"""
        if not self.__od_amf_orig:
            # if not self.optical_depth_amf_offsets['460.3']:
            #     txt = 'please define an od offset (miniSASP only measures relative od) by setting optical_depth_amf_offset'
            #     raise AttributeError(txt)
            self.__od_amf_orig = timeseries.TimeSeries(-1 * np.log(self.sun_intensities.data))
            self.__od_amf_orig._data_period = self.sun_intensities._data_period

        if self.optical_depth_amf_offsets != self.__od_afm_offset_last:
            self.__od_amf = self.__od_amf_orig
            cols = self.__od_amf.data.columns
            for e, col in enumerate(cols):
                self.__od_amf.data[col] += self.optical_depth_amf_offsets[col]
        return self.__od_amf

    @property
    def optical_depth(self):
        """optical depth + unkonwn offset. after determining the offset you might want to set the property sup_offsets"""
        if not self.__od:
            if not self.optical_depth_amf_offsets['460.3']:
                txt = 'please define an od offset (miniSASP only measures relative od) by setting optical_depth_amf_offset'
                raise AttributeError(txt)
            self.__od = self.optical_depth_amf / self.air_mass_factor
        return self.__od


    @property
    def optical_depth_rayleigh(self):
        if not self.__housekeeping:
            txt = 'To claculate the optical depth you will have to set the housekeeping property with a timeseries containing Altitude (m) Temperature (C) and Pressure (hPa)'
            raise AttributeError(txt)

        if not self.__od_ray_orig:
            hkt = self.housekeeping.copy()
            # import pdb
            # pdb.set_trace()
            # to make sure there is no two values with the same altitude in a row (causes sims to fail) ...
            clean = False
            while not clean:
                alt = hkt.data.Altitude.values
                dist = alt[1:] - alt[:-1]
                if (dist == 0).sum():
                    hkt.data.iloc[np.where(dist == 0)] = np.nan
                else:
                    clean = True
            # pdb.set_trace()
            hkt.data.dropna(axis=1, inplace=True)
            # pdb.set_trace()
            # do the od calculations for each wavelength channel
            c_list = np.array(sorted(self.wavelength_channels.items(), key=lambda x: x[1]))[:, 0]
            for e, c in enumerate(c_list):
                wl = self.wavelength_channels[c]

                out = []
                for i in range(hkt.data.shape[0]):
                    if i < 1:
                        out.append(0.)
                        continue
                    hktt = hkt.data.iloc[0:i, :]
                    alt = -1 * hktt.Altitude.values
                    press = hktt.Pressure_Pa.values
                    temp = hktt.Temperature.values
                    ray = _bray.rayleigh_optical_depth(alt, press, temp + 273.15, float(wl))
                    out.append(ray)
                # pdb.set_trace()
                hkt.data[wl] = np.array(out)
                # pdb.set_trace()
            wll = np.array(sorted(self.wavelength_channels.items(), key=lambda x: x[1]))[:, 1]#.astype(float)
            # pdb.set_trace()
            hkt = hkt._del_all_columns_but(wll)
            # pdb.set_trace()
            self.__od_ray_orig = hkt


        if self.optical_depth_rayleigh_offsets != self.__od_ray_offset_last:
            # pdb.set_trace()
            self.__od_ray = self.__od_ray_orig.copy()
            # pdb.set_trace()
            # apply offset
            wll = np.array(sorted(self.wavelength_channels.items(), key=lambda x: x[1]))[:, 1]#.astype(float)
            for wl in wll:
                self.__od_ray.data[wl] += self.optical_depth_rayleigh_offsets[wl]

            self.__od_ray_offset_last = self.optical_depth_rayleigh_offsets
        # pdb.set_trace()
        return self.__od_ray

    @optical_depth_rayleigh.setter
    def optical_depth_rayleigh(self, ts):
        self.__od_ray = ts
        self.__od_ray_orig = ts

    @property
    def aerosol_optical_depth(self):
        if not self.__aod:
            self.__aod = self.optical_depth - self.optical_depth_rayleigh
        return self.__aod






class Sun_Intensities_TS(timeseries.TimeSeries):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self._data_period = 30
    # def plot(self, offset=[0, 0, 0, 0], airmassfct=True, move_max=True, legend=True, all_on_one_axis = False,
    #          additional_axes=False,
    #          errors = False,
    #          rayleigh=True):
    #     """plots ... sorry, but this is a messi function. Things should have been done different, e.g too much data
    #      processing whith the data not put out ... need fixn
    #     Arguments
    #     ---------
    #     offset: list
    #     airmassfct: bool.
    #         If the airmass factor is included or not.
    #         True: naturally the air-mass factor is included in the data, so this does nothing.
    #         False: data is corrected to correct for the slant angle
    #     rayleigh: bool or the aod part of the output of miniSASP.simulate_from_size_dist_LS.
    #         make sure there is no airmassfkt included in this!!
    #     all_on_one_axis: bool or axes instance
    #         if True all is plotted in one axes. If axes instances this axis is used.
    #     """
    #
    #     m_size = 5
    #     m_ewidht = 1.5
    #     l_width = 2
    #     gridspec_kw = {'wspace': 0.05}
    #     no_axes = 4
    #     if all_on_one_axis:
    #         no_axes = 1
    #     if additional_axes:
    #         no_axes = no_axes + additional_axes
    #
    #     if type(all_on_one_axis).__name__ == 'AxesSubplot':
    #         a = all_on_one_axis
    #         f = a.get_figure()
    #     else:
    #         f, a = plt.subplots(1, no_axes, gridspec_kw=gridspec_kw)
    #     columns = ['460.3', '460.3 max', '550.4', '550.4 max', '671.2', '671.2 max', '860.7', '860.7 max']
    #     # peaks_max = [460.3, '460.3 max', 550.4, '550.4 max', 860.7, '860.7 max', 671.2,
    #     #        '671.2 max']
    #     if not all_on_one_axis:
    #         f.set_figwidth(15)
    #     #################
    #     for i in range(int(len(columns) / 2)):
    #         col = plt_tools.wavelength_to_rgb(columns[i * 2]) * 0.8
    #         intens = self.data[columns[i * 2]].dropna()  # .plot(ax = a, style = 'o', label = '%s nm'%colums[i*2])
    #         x = intens.index.get_level_values(1)
    #         if type(rayleigh) == bool:
    #             if rayleigh:
    #                 rayleigh_corr = 0
    #         else:
    #             # print('mach ick')
    #             aodt = rayleigh[float(columns[i * 2])].loc[:, ['rayleigh']]
    #             intenst = intens.copy()
    #             intenst.index = intenst.index.droplevel(['Time', 'Sunelevation'])
    #             aodt_sit = pd.concat([aodt, intenst]).sort_index().interpolate()
    #             aodt_sit = aodt_sit.groupby(aodt_sit.index).mean().reindex(intenst.index)
    #             rayleigh_corr = aodt_sit.rayleigh.values / np.sin(intens.index.get_level_values(2))
    #             # return aodt
    #
    #         if not airmassfct:
    #             amf_corr = np.sin(intens.index.get_level_values(2))
    #         else:
    #             amf_corr = 1
    #         if not all_on_one_axis:
    #             atmp = a[i]
    #         else:
    #             atmp = a
    #
    #
    #         y = (offset[i] - np.log(intens) - rayleigh_corr) * amf_corr
    #         g, = atmp.plot(y, x)
    #         g.set_label('%s nm' % columns[i * 2])
    #         g.set_linestyle('')
    #         g.set_marker('o')
    #         #         g = a.get_lines()[-1]
    #         g.set_markersize(m_size)
    #         g.set_markeredgewidth(m_ewidht)
    #         g.set_markerfacecolor('None')
    #         g.set_markeredgecolor(col)
    #
    #         if move_max:
    #             #             sun_intensities.data.iloc[:,i*2+1].dropna().plot(ax = a)
    #             intens = self.data[
    #                 columns[i * 2 + 1]].dropna()  # .plot(ax = a, style = 'o', label = '%s nm'%colums[i*2])
    #             x = intens.index.values
    #
    #             g, = a[i].plot(offset[i] - np.log(intens), x)
    #             #             g = a.get_lines()[-1]
    #             g.set_color(col)
    #             #             g.set_solid_joinstyle('round')
    #             g.set_linewidth(l_width)
    #             g.set_label(None)
    #
    #         if i != 0 and not all_on_one_axis:
    #             atmp.set_yticklabels([])
    #
    #         if i == 4:
    #             break
    #     if all_on_one_axis:
    #         a.legend()
    #     else:
    #         if legend:
    #             for aa in a:
    #                 aa.legend()
    #     if not airmassfct:
    #         txt = 'OD'
    #     else:
    #         txt = 'OD * (air-mass factor)'`
    #     if all_on_one_axis:
    #         atmp = a
    #     else:
    #         atmp = a[0]
    #     atmp.set_xlabel(txt)
    #     if not all_on_one_axis:
    #         atmp.xaxis.set_label_coords(2.05, -0.07)
    #     atmp.set_ylabel('Altitude (m)')
    #     return a


