import csv
import io
from pathlib import Path
import struct
import tempfile
import unittest

import numpy as np

from atmPy.io.campbell.tob3 import TOB3Reader, read_tob3


LOGGER_FILE = Path(
    "/Users/htelg/data/arm/raw_from_SDcard/2m/202608_afterLightningStrike/"
    "49238.SEARM_RADSYS_1Hz.dat"
)


def _csv_line(values):
    output = io.StringIO(newline="")
    csv.writer(output, lineterminator="\r\n", quoting=csv.QUOTE_ALL).writerow(values)
    return output.getvalue().encode("ascii")


def _write_small_tob3(path):
    validation_stamp = 12_345
    # One padding byte follows each 8-byte record.  Campbell frames may divide
    # otherwise unused body bytes evenly between records.
    frame_size = 34
    header = b"".join(
        [
            _csv_line(
                ["TOB3", "station", "CR1000X", "123", "OS", "CPU:test", "1", "created"]
            ),
            _csv_line(["table", "1 SEC", frame_size, 4, validation_stamp, "Sec100Usec"]),
            _csv_line(["temperature", "flag", "humidity"]),
            _csv_line(["DegC", "", "%"]),
            _csv_line(["Smp", "Smp", "Smp"]),
            _csv_line(["IEEE4B", "UINT2", "FP2"]),
        ]
    )
    frames = []
    for frame_index in range(2):
        first_second = 1 + frame_index * 2
        first_record = 10 + frame_index * 2
        body = struct.pack("<III", first_second, 0, first_record)
        for in_frame in range(2):
            temperature = 20.0 + frame_index * 2 + in_frame
            flag = frame_index * 2 + in_frame
            fp2_humidity = 0x4000 | (5000 + frame_index * 10 + in_frame)
            body += struct.pack(">fHHx", temperature, flag, fp2_humidity)
        body += struct.pack("<HH", 0, validation_stamp)
        frames.append(body)
    path.write_bytes(header + b"".join(frames))


class TOB3ReaderTest(unittest.TestCase):
    def test_vectorized_eager_decode_and_subsetting(self):
        with tempfile.TemporaryDirectory() as directory:
            path = Path(directory) / "small.dat"
            _write_small_tob3(path)

            ds = read_tob3(
                path,
                time=slice("1990-01-01T00:00:02", "1990-01-01T00:00:04"),
                variables=["temperature", "flag", "humidity"],
                chunks=None,
            )
            output_path = Path(directory) / "subset.nc"
            ds.to_netcdf(output_path)
            self.assertGreater(output_path.stat().st_size, 0)

        np.testing.assert_array_equal(
            ds.datetime.data,
            np.array(["1990-01-01T00:00:02", "1990-01-01T00:00:03"], dtype="datetime64[ns]"),
        )
        np.testing.assert_array_equal(ds.record_number.data, [11, 12])
        np.testing.assert_allclose(ds.temperature.data, [21.0, 22.0])
        np.testing.assert_array_equal(ds.flag.data, [1, 2])
        np.testing.assert_allclose(ds.humidity.data, [50.01, 50.10], rtol=1e-6)
        self.assertEqual(ds.temperature.dims, ("datetime",))
        self.assertIn("datetime", ds.indexes)
        self.assertEqual(ds.temperature.attrs["units"], "DegC")
        plot = ds.temperature.plot()
        self.assertEqual(plot[0].axes.get_xlabel(), "datetime")
        self.assertTrue(np.issubdtype(plot[0].get_xdata().dtype, np.datetime64))

    @unittest.skipUnless(LOGGER_FILE.exists(), "2 GB logger regression file is not available")
    def test_large_logger_file_header_ring_order_and_values(self):
        reader = TOB3Reader(LOGGER_FILE)

        self.assertEqual(reader.header.table_name, "SEARM_RADSYS_1Hz")
        self.assertEqual(reader.header.frame_size, 1006)
        self.assertEqual(reader.records_per_frame, 9)
        self.assertEqual(reader.record_count, 19_207_732)
        self.assertEqual(reader.start_time, np.datetime64("2026-01-05T07:08:08", "ns"))
        self.assertEqual(reader.end_time, np.datetime64("2026-08-18T14:50:10", "ns"))

        first = reader.read(
            records=slice(0, 2),
            variables=["LoggerPS_Bat", "CleanFlag", "SWd", "Ta", "RH"],
            chunks=None,
        )
        np.testing.assert_array_equal(first.record_number.data, [17_308_872, 17_308_873])
        np.testing.assert_allclose(
            first.LoggerPS_Bat.data,
            [14.203951835632324, 14.205350875854492],
        )
        np.testing.assert_array_equal(first.CleanFlag.data, [0, 0])
        np.testing.assert_allclose(first.SWd.data, [0.0048440509, 0.0048668939])
        np.testing.assert_allclose(first.Ta.data, [3.464, 3.491], rtol=1e-6)
        np.testing.assert_allclose(first.RH.data, [83.5, 83.5], rtol=1e-6)

    @unittest.skipUnless(LOGGER_FILE.exists(), "2 GB logger regression file is not available")
    def test_large_logger_file_lazy_time_subset_stops_at_minor_footer(self):
        ds = read_tob3(
            LOGGER_FILE,
            time=slice("2026-08-18T14:50:07", "2026-08-18T14:50:12"),
            variables=["Ta", "RH"],
            chunks=2,
        )

        self.assertTrue(hasattr(ds.Ta.data, "chunks"))
        self.assertTrue(hasattr(ds.datetime.data, "chunks"))
        self.assertEqual(ds.Ta.dims, ("datetime",))
        self.assertNotIn("datetime", ds.indexes)
        loaded = ds.compute()
        np.testing.assert_array_equal(
            loaded.datetime.data,
            np.array(
                [
                    "2026-08-18T14:50:07",
                    "2026-08-18T14:50:08",
                    "2026-08-18T14:50:09",
                    "2026-08-18T14:50:10",
                ],
                dtype="datetime64[ns]",
            ),
        )
        np.testing.assert_allclose(loaded.Ta.data, [24.38, 24.38, 24.38, 24.39], rtol=1e-6)
        np.testing.assert_allclose(loaded.RH.data, [94.0, 94.1, 94.1, 94.1], rtol=1e-6)


if __name__ == "__main__":
    unittest.main()
