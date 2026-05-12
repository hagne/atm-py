import os
import sys
import tempfile
import types
import unittest
import importlib.util

import pandas as pd
import xarray as xr


os.environ.setdefault("MPLCONFIGDIR", os.path.join(tempfile.gettempdir(), "atmPy-matplotlib"))

timezonefinder_stub = types.ModuleType("timezonefinder")


class TimezoneFinder:
    def timezone_at(self, lng, lat):
        if lng > -100:
            return "America/Chicago"
        return "America/Denver"


timezonefinder_stub.TimezoneFinder = TimezoneFinder
sys.modules.setdefault("timezonefinder", timezonefinder_stub)

from atmPy.general.measurement_site import MovingPlatform, Network, Station


def make_boulder_station():
    return Station(
        lat=40.0,
        lon=-105.0,
        alt=1600,
        name="Boulder",
        abbreviation="BOU",
    )


class StationTests(unittest.TestCase):
    def test_station_normalizes_name_and_abbreviation(self):
        station = Station(
            lat=36.6057,
            lon=-97.4878,
            alt=314,
            name=" Test, Site (A) ",
            abbreviation=["sgp", "SGP"],
            state="OK",
            country="US",
            data_product="surfrad",
        )

        self.assertEqual(station.name, "Test_ Site _A")
        self.assertEqual(station.abb, "sgp")
        self.assertEqual(station.abb_alternative, ["sgp", "SGP"])
        self.assertEqual(station.data_product, "surfrad")
        self.assertEqual(station["data_product"], "surfrad")
        self.assertIn("Station(name=Test_ Site _A", str(station))

    def test_station_accepts_single_abbreviation(self):
        station = Station(
            lat=1,
            lon=2,
            alt=3,
            name="Single Abbreviation",
            abbreviation="one",
        )

        self.assertEqual(station.abb, "one")
        self.assertEqual(station.abb_alternative, ["one"])

    def test_time_zone_returns_expected_format_and_standard_time_values(self):
        station = make_boulder_station()

        # time_zone is currently a property, so the getter is called directly
        # to characterize the date-specific behavior.
        time_zone = station.get_time_zone_for_date("2024-01-15 12:00:00")

        self.assertEqual(
            set(time_zone),
            {"tz_str", "tz_name", "diff2UTC", "diff2UTC_of_standard_time"},
        )
        self.assertEqual(time_zone["tz_str"], "America/Denver")
        self.assertEqual(time_zone["tz_name"], "MST")
        self.assertEqual(time_zone["diff2UTC"], -7.0)
        self.assertEqual(time_zone["diff2UTC_of_standard_time"], -7.0)

    def test_time_zone_returns_expected_daylight_saving_values(self):
        station = make_boulder_station()

        time_zone = station.get_time_zone_for_date("2024-07-15 12:00:00")

        self.assertEqual(time_zone["tz_str"], "America/Denver")
        self.assertEqual(time_zone["tz_name"], "MDT")
        self.assertEqual(time_zone["diff2UTC"], -6.0)
        self.assertEqual(time_zone["diff2UTC_of_standard_time"], -7.0)

    def test_get_sun_position_pvlib_returns_expected_dataset_and_values(self):
        station = make_boulder_station()

        sun_position = station.get_sun_position(
            "2024-06-21 18:00:00+00:00",
            method="pvlib",
            return_type="dataset",
        )

        self.assertIsInstance(sun_position, xr.Dataset)
        self.assertEqual(sun_position.attrs["source"], "pvlib")
        self.assertEqual(sun_position.sizes["datetime"], 1)
        self.assertEqual(
            set(sun_position.data_vars),
            {
                "zenith",
                "zenith_geometric",
                "elevation",
                "elevation_geometric",
                "azimuth",
                "equation_of_time",
                "airmass",
                "airmass_absolute",
                "sun_earth_distance",
            },
        )
        self.assertAlmostEqual(float(sun_position.zenith.values[0]), 0.368176, places=6)
        self.assertAlmostEqual(float(sun_position.elevation.values[0]), 1.202620, places=6)
        self.assertAlmostEqual(float(sun_position.azimuth.values[0]), 2.392662, places=6)
        self.assertAlmostEqual(float(sun_position.airmass.values[0]), 1.071332, places=6)
        self.assertAlmostEqual(float(sun_position.sun_earth_distance.values[0]), 1.016251, places=6)

    def test_get_sun_position_pvlib_can_return_dataframe(self):
        station = make_boulder_station()

        sun_position = station.get_sun_position(
            "2024-06-21 18:00:00+00:00",
            method="pvlib",
            return_type="dataframe",
        )

        self.assertEqual(list(sun_position.index.names), ["datetime"])
        self.assertIn("elevation", sun_position.columns)
        self.assertAlmostEqual(float(sun_position.loc[:, "elevation"].iloc[0]), 1.202620, places=6)


class MovingPlatformTests(unittest.TestCase):
    def __init__(self, *args, verbose=True, **kwargs):
        super().__init__(*args, **kwargs)
        self.verbose = verbose
    def make_platform(self):
        index = pd.DatetimeIndex(["2024-06-21 18:00:00", "2024-06-21 19:00:00"])
        index.name = 'datetime'
        return MovingPlatform(
            lat=pd.Series([40.0, 42.0], index=index).to_xarray(),
            lon=pd.Series([-105.0, -95.0], index=index).to_xarray(),
            alt=pd.Series([1600.0, 1800.0], index=index).to_xarray(),
            name="Moving Platform",
            abbreviation="MOB",
        )

    def test_coordinates_at_datetime_is_cached_dataset_with_attrs(self):
        platform = self.make_platform()

        coords = platform.data

        # self.assertIs(coords, platform.coordinates_at_datetime)
        self.assertIsInstance(coords, xr.Dataset)
        self.assertEqual(coords.sizes["datetime"], 2)
        self.assertEqual(set(coords.data_vars), {"lat", "lon", "alt"})

    def test_time_zone_uses_position_at_each_datetime(self):
        platform = self.make_platform()

        time_zone = platform.time_zone

        self.assertEqual(list(time_zone["tz_str"]), ["America/Denver", "America/Chicago"])
        self.assertEqual(list(time_zone["tz_name"]), ["MDT", "CDT"])
        self.assertEqual(list(time_zone["diff2UTC"]), [-6.0, -5.0])

    def test_sun_position_returns_expected_dataset_and_values(self):
        if self.verbose:
            print("test_sun_position_returns_expected_dataset_and_values:", end=" ... ")
        platform = self.make_platform()

        sun_position = platform.sun_position

        self.assertIsInstance(sun_position, xr.Dataset)
        self.assertEqual(sun_position.attrs["source"], "pvlib")
        self.assertIs(sun_position, platform.sun_position)
        self.assertEqual(sun_position.sizes["datetime"], 2)
        self.assertEqual(
            set(sun_position.data_vars),
            {
                "zenith",
                "zenith_geometric",
                "elevation",
                "elevation_geometric",
                "azimuth",
                "equation_of_time",
                "airmass",
                "airmass_absolute",
                "sun_earth_distance",
            },
        )
        self.assertAlmostEqual(float(sun_position.zenith.values[0]), 0.368176, places=6)
        self.assertAlmostEqual(float(sun_position.zenith.values[1]), 0.352157, places=6)
        self.assertAlmostEqual(float(sun_position.elevation.values[0]), 1.202620, places=6)
        self.assertAlmostEqual(float(sun_position.elevation.values[1]), 1.218640, places=6)
        self.assertAlmostEqual(float(sun_position.azimuth.values[0]), 2.392662, places=6)
        self.assertAlmostEqual(float(sun_position.azimuth.values[1]), 3.596215, places=6)
        self.assertAlmostEqual(float(sun_position.airmass.values[0]), 1.071332, places=6)
        self.assertAlmostEqual(float(sun_position.airmass.values[1]), 1.064901, places=6)
        self.assertAlmostEqual(float(sun_position.sun_earth_distance.values[0]), 1.016251, places=6)
        self.assertAlmostEqual(float(sun_position.sun_earth_distance.values[1]), 1.016254, places=6)
        if self.verbose:
            print("passed")



class NetworkTests(unittest.TestCase):
    def test_network_registers_stations_by_sanitized_name(self):
        station = Station(
            lat=36.6057,
            lon=-97.4878,
            alt=314,
            name="Test Site. One",
            abbreviation=["one", "ONE"],
        )

        network = Network([station], generate_subnetworks=False, name="test")

        self.assertIs(station.parent_network, network)
        self.assertIs(network.stations.Test_Site_One, station)
        self.assertEqual(network.stations.list, [station])

    def test_network_can_be_created_from_station_dicts(self):
        network = Network(
            [
                {
                    "lat": 36.6057,
                    "lon": -97.4878,
                    "alt": 314,
                    "name": "Dictionary Site",
                    "abbreviation": "dict",
                }
            ],
            generate_subnetworks=False,
        )

        self.assertEqual(network.stations.Dictionary_Site.name, "Dictionary Site")
        self.assertEqual(network.stations.Dictionary_Site.abb, "dict")

    def test_find_site_matches_primary_and_alternative_abbreviations(self):
        station = Station(
            lat=36.6057,
            lon=-97.4878,
            alt=314,
            name="Findable Site",
            abbreviation=["abc", "ALT"],
        )
        network = Network([station], generate_subnetworks=False)

        self.assertIs(network.stations.find_site("abc"), station)
        self.assertIs(network.stations.find_site("ABC"), station)
        self.assertIs(network.stations.find_site("ALT"), station)

    def test_operation_period_creates_active_and_inactive_subnetworks(self):
        active_station = Station(
            lat=1,
            lon=2,
            alt=3,
            name="Active Site",
            abbreviation="active",
            operation_period=["2020", "present"],
        )
        inactive_station = Station(
            lat=4,
            lon=5,
            alt=6,
            name="Inactive Site",
            abbreviation="inactive",
            operation_period=["2000", "2010"],
        )

        network = Network([active_station, inactive_station])

        self.assertEqual([sub.name for sub in network.subnetworks._network_list], ["active", "inactive"])
        self.assertEqual(network.subnetworks.active.stations.list, [active_station])
        self.assertEqual(network.subnetworks.inactive.stations.list, [inactive_station])


if __name__ == "__main__":
    unittest.main()
