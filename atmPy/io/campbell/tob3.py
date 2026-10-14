"""Fast, subset-aware reader for Campbell Scientific TOB3 files.

TOB3 is the binary, fixed-frame format written to datalogger storage cards.
Each file starts with six CSV header lines.  The remaining bytes are frames
containing a 12-byte timestamp/record header, several data records, and a
4-byte validation footer.  Storage cards are circular: after the allocated
space fills, new frames overwrite the physical beginning of the file.

The implementation follows Campbell Scientific's ``TOB2 or TOB3`` and
``Binary Data Value Types`` documentation in the LoggerNet manual.  It avoids
creating one Python object per value: complete frame ranges are decoded with
structured NumPy dtypes, while dask can keep large reads lazy and bounded.

Typical use
-----------
Read a small record range eagerly::

    ds = read_tob3(path, records=slice(0, 1_000), chunks=None)
    ds.to_netcdf("first_1000_records.nc")

Read a time range and only selected fields::

    ds = read_tob3(
        path,
        time=slice("2026-08-18T14:00", "2026-08-18T15:00"),
        variables=["SWd", "Ta", "RH"],
        chunks=100_000,
    )

``chunks`` is a number of records per dask chunk.  Its non-None default makes
opening a multi-gigabyte file safe; use ``chunks=None`` for an eager subset.
The time slice is half-open: its start is included and its stop is excluded.

References
----------
Campbell Scientific, `TOB2 or TOB3
<https://help.campbellsci.com/loggernet-manual/ln_manual/
campbell_scientific_file_formats/tob2_or_tob3.htm>`_ and `Binary Data Value
Types <https://help.campbellsci.com/loggernet-manual/ln_manual/
campbell_scientific_file_formats/binary_data_value_types.htm>`_.
"""

from __future__ import annotations

import csv
import dataclasses
import datetime as dt
import decimal
import operator
from collections.abc import Iterable, Iterator
from pathlib import Path
import re
import struct

import numpy as np
import xarray as xr

_FRAME_HEADER_SIZE = 12
_FRAME_FOOTER_SIZE = 4
_CAMPBELL_EPOCH_NS = np.datetime64("1990-01-01T00:00:00", "ns").astype(np.int64)
_FRAME_HEADER = struct.Struct("<III")
_FRAME_FOOTER = struct.Struct("<HH")
_ASCII_TYPE = re.compile(r"ASCII\((\d+)\)", re.IGNORECASE)
_FRAME_RESOLUTION = re.compile(r"Sec(\d*)(Nsec|Usec|Msec)", re.IGNORECASE)


@dataclasses.dataclass(frozen=True)
class TOB3Field:
    """Description of one measured field from TOB3 header lines 3--6."""

    name: str
    units: str
    processing: str
    data_type: str


@dataclasses.dataclass(frozen=True)
class TOB3Header:
    """Parsed metadata and binary layout from the six TOB3 header lines."""

    file_type: str
    station_name: str
    logger_model: str
    logger_serial_number: str
    operating_system: str
    program_name: str
    program_signature: str
    file_created: str
    table_name: str
    record_interval: str
    frame_size: int
    intended_table_size: int
    validation_stamp: int
    frame_time_resolution: str
    extra_table_metadata: tuple[str, ...]
    fields: tuple[TOB3Field, ...]
    data_offset: int


@dataclasses.dataclass(frozen=True)
class _FrameInfo:
    seconds: int
    subsecond_ticks: int
    first_record_number: int
    footer: int
    footer_validation: int


def read_tob3_header(path: str | Path) -> TOB3Header:
    """Read and validate the six ASCII header lines of a TOB3 file.

    Only the small header is read.  This is useful for inspecting an unknown
    file before deciding which variables or record interval to load.

    Parameters
    ----------
    path
        Campbell Scientific TOB3 file.

    Returns
    -------
    TOB3Header
        Immutable, typed representation of the file metadata.
    """
    path = Path(path)
    rows: list[list[str]] = []
    with path.open("rb") as file:
        for line_number in range(1, 7):
            line = file.readline()
            if not line:
                raise ValueError(f"{path} ends before TOB3 header line {line_number}")
            try:
                text = line.decode("ascii")
            except UnicodeDecodeError as error:
                raise ValueError(
                    f"TOB3 header line {line_number} in {path} is not ASCII"
                ) from error
            rows.append([value.strip() for value in next(csv.reader([text]))])
        data_offset = file.tell()

    if len(rows[0]) < 8 or rows[0][0].upper() != "TOB3":
        file_type = rows[0][0] if rows[0] else ""
        raise ValueError(f"Expected a TOB3 file, found {file_type!r} in {path}")
    if len(rows[1]) < 6:
        raise ValueError(f"TOB3 table header in {path} has fewer than six fields")

    names, units, processing, data_types = rows[2:6]
    field_counts = {len(names), len(units), len(processing), len(data_types)}
    if len(field_counts) != 1:
        raise ValueError(
            "TOB3 field names, units, processing, and data types have different lengths"
        )
    if len(set(names)) != len(names):
        raise ValueError("TOB3 field names must be unique")

    fields = tuple(
        TOB3Field(name, unit, process, data_type)
        for name, unit, process, data_type in zip(
            names, units, processing, data_types, strict=True
        )
    )
    environment = rows[0]
    table = rows[1]
    return TOB3Header(
        file_type=environment[0],
        station_name=environment[1],
        logger_model=environment[2],
        logger_serial_number=environment[3],
        operating_system=environment[4],
        program_name=environment[5],
        program_signature=environment[6],
        file_created=environment[7],
        table_name=table[0],
        record_interval=table[1],
        frame_size=int(table[2]),
        intended_table_size=int(table[3]),
        validation_stamp=int(table[4]),
        frame_time_resolution=table[5],
        extra_table_metadata=tuple(table[6:]),
        fields=fields,
        data_offset=data_offset,
    )


class TOB3Reader:
    """Inspect and read one Campbell Scientific TOB3 file.

    Construction reads only the six-line header and a small number of frame
    headers/footers.  In particular, it does not scan or load the whole file.
    The physical-to-chronological frame rotation is found with a binary search,
    which makes opening a full circular file fast even when it is several GB.

    Notes
    -----
    TOB3 timestamps do not carry a time-zone offset.  Returned ``datetime64``
    values therefore represent the logger's naive clock exactly as stored.
    """

    def __init__(self, path: str | Path):
        self.path = Path(path)
        self.header = read_tob3_header(self.path)
        self._record_dtype = _record_dtype(self.header.fields)
        self.record_size = self._record_dtype.itemsize
        self.records_per_frame = (
            self.header.frame_size - _FRAME_HEADER_SIZE - _FRAME_FOOTER_SIZE
        ) // self.record_size
        if self.records_per_frame < 1:
            raise ValueError(
                f"Frame size {self.header.frame_size} cannot contain a "
                f"{self.record_size}-byte data record"
            )

        frame_body_size = self.header.frame_size - _FRAME_HEADER_SIZE - _FRAME_FOOTER_SIZE
        padding_bytes = frame_body_size - self.records_per_frame * self.record_size
        padding_per_record, padding_remainder = divmod(
            padding_bytes, self.records_per_frame
        )
        if padding_remainder:
            raise ValueError(
                f"TOB3 frame padding ({padding_bytes} bytes) cannot be divided "
                f"between {self.records_per_frame} records"
            )
        self.record_stride = self.record_size + padding_per_record
        self._frame_record_dtype = _dtype_with_itemsize(
            self._record_dtype, self.record_stride
        )

        payload_size = self.path.stat().st_size - self.header.data_offset
        self.physical_frame_count, trailing_bytes = divmod(
            payload_size, self.header.frame_size
        )
        if trailing_bytes:
            raise ValueError(
                f"{self.path} has {trailing_bytes} bytes after its last complete frame"
            )
        if self.physical_frame_count == 0:
            raise ValueError(f"{self.path} contains no TOB3 data frames")

        self._frame_timestamp_tick_ns = _parse_frame_resolution(
            self.header.frame_time_resolution
        )
        self._record_interval_ns = _parse_record_interval(self.header.record_interval)
        self.frame_count = self._count_valid_frames()
        self._rotation = self._find_rotation()
        self._last_frame_record_count = self._read_last_frame_record_count()
        self.record_count = (
            (self.frame_count - 1) * self.records_per_frame
            + self._last_frame_record_count
        )

    @property
    def start_time(self) -> np.datetime64:
        """Timestamp of the oldest retained record."""
        return np.datetime64(self._timestamp_ns_at(0), "ns")

    @property
    def end_time(self) -> np.datetime64:
        """Timestamp of the newest retained record."""
        return np.datetime64(self._timestamp_ns_at(self.record_count - 1), "ns")

    def read(
        self,
        *,
        records: slice | None = None,
        time: slice | None = None,
        variables: Iterable[str] | None = None,
        chunks: int | None = 250_000,
    ) -> xr.Dataset:
        """Read a chronological subset as an :class:`xarray.Dataset`.

        Parameters
        ----------
        records
            Zero-based chronological record-position slice.  Python negative
            bounds are supported; the step must be one.  This is distinct from
            the datalogger's ``record_number`` coordinate.
        time
            Half-open datetime slice ``[start, stop)``.  Strings, Python
            datetimes, and NumPy datetimes are accepted.  The lookup uses frame
            timestamps and does not assume the file has no acquisition gaps.
        variables
            Exact field names to decode.  By default all fields are returned.
        chunks
            Number of records per lazy dask chunk.  The default keeps large
            files out of memory.  Pass ``None`` to eagerly read a small subset.

        Returns
        -------
        xarray.Dataset
            Variables use ``datetime`` as their dimension, so xarray plots
            against time by default.  ``record_number`` is an auxiliary
            coordinate.  For chunked reads the datetime coordinate remains
            lazy and is not turned into an in-memory xarray index.

        Raises
        ------
        ValueError
            If both record and time slices are supplied, a requested variable
            is absent, or the selected frames fail validation.
        """
        start, stop = self._normalize_selection(records=records, time=time)
        selected_fields = self._select_fields(variables)

        lazy = chunks is not None and start != stop
        if not lazy:
            arrays = self._read_record_range(start, stop, selected_fields)
        else:
            if not isinstance(chunks, int) or isinstance(chunks, bool) or chunks <= 0:
                raise ValueError("chunks must be a positive record count or None")
            arrays = self._lazy_record_range(start, stop, selected_fields, chunks)

        data_vars = {
            field.name: (
                "datetime",
                arrays[field.name],
                {
                    "units": field.units,
                    "processing": field.processing,
                    "campbell_data_type": field.data_type,
                },
            )
            for field in selected_fields
        }
        coordinate_values = {
            "datetime": (
                "datetime",
                arrays["datetime"],
                {"timezone": "logger local time; offset not stored by TOB3"},
            ),
            "record_number": (
                "datetime",
                arrays["record_number"],
                {"long_name": "Campbell datalogger record number"},
            ),
        }
        if lazy:
            # Creating the usual pandas-backed xarray index would compute the
            # entire dask datetime coordinate while opening the dataset.
            coordinates = xr.Coordinates(coordinate_values, indexes={})
        else:
            coordinates = coordinate_values
        return xr.Dataset(
            data_vars=data_vars,
            coords=coordinates,
            attrs=self._dataset_attrs(),
        )

    def iter_chunks(
        self,
        *,
        records: slice | None = None,
        time: slice | None = None,
        variables: Iterable[str] | None = None,
        records_per_chunk: int = 250_000,
    ) -> Iterator[xr.Dataset]:
        """Yield eager, chronological datasets with bounded memory use.

        This is useful for streaming processing or writing one subset at a
        time.  Selection semantics match :meth:`read`.
        """
        if (
            not isinstance(records_per_chunk, int)
            or isinstance(records_per_chunk, bool)
            or records_per_chunk <= 0
        ):
            raise ValueError("records_per_chunk must be a positive integer")
        start, stop = self._normalize_selection(records=records, time=time)
        for chunk_start in range(start, stop, records_per_chunk):
            chunk_stop = min(chunk_start + records_per_chunk, stop)
            yield self.read(
                records=slice(chunk_start, chunk_stop),
                variables=variables,
                chunks=None,
            )

    def _count_valid_frames(self) -> int:
        if not self._frame_is_valid(0):
            raise ValueError("The first TOB3 frame does not have a valid footer stamp")
        if self._frame_is_valid(self.physical_frame_count - 1):
            return self.physical_frame_count

        # A new card is preallocated with invalid frames.  Until the first
        # wrap, the written frames form one valid prefix, found in O(log n).
        low, high = 1, self.physical_frame_count
        while low < high:
            middle = (low + high) // 2
            if self._frame_is_valid(middle):
                low = middle + 1
            else:
                high = middle
        return low

    def _find_rotation(self) -> int:
        """Return the physical index of the oldest valid frame."""
        if self.frame_count == 1:
            return 0
        first = self._frame_sort_key(0)
        last = self._frame_sort_key(self.frame_count - 1)
        if first <= last:
            return 0

        # Standard minimum search in an ascending array rotated once by the
        # datalogger's circular overwrite operation.
        low, high = 0, self.frame_count - 1
        while low < high:
            middle = (low + high) // 2
            if self._frame_sort_key(middle) > self._frame_sort_key(high):
                low = middle + 1
            else:
                high = middle
        return low

    def _read_last_frame_record_count(self) -> int:
        physical_index = (self._rotation + self.frame_count - 1) % self.frame_count
        info = self._read_frame_info(physical_index)
        footer_offset = info.footer & 0x07FF
        if footer_offset == 0:
            return self.records_per_frame

        # In a minor (partially written) frame, the footer at the physical end
        # stores the backward byte offset to the embedded footer immediately
        # following the last complete data record.
        embedded_footer_position = (
            self.header.frame_size - _FRAME_FOOTER_SIZE - footer_offset
        )
        data_bytes = embedded_footer_position - _FRAME_HEADER_SIZE
        record_count, remainder = divmod(data_bytes, self.record_stride)
        if remainder or not 0 <= record_count <= self.records_per_frame:
            raise ValueError("The final TOB3 minor-frame offset is inconsistent")

        with self.path.open("rb") as file:
            file.seek(
                self.header.data_offset
                + physical_index * self.header.frame_size
                + embedded_footer_position
            )
            footer_bytes = file.read(_FRAME_FOOTER_SIZE)
        if len(footer_bytes) != _FRAME_FOOTER_SIZE:
            raise ValueError("The final TOB3 minor-frame footer is truncated")
        _, validation = _FRAME_FOOTER.unpack(footer_bytes)
        if validation not in self._valid_footer_stamps:
            raise ValueError("The final TOB3 minor-frame footer is invalid")
        return record_count

    @property
    def _valid_footer_stamps(self) -> tuple[int, int]:
        stamp = self.header.validation_stamp
        return stamp, 0xFFFF - stamp

    def _frame_is_valid(self, physical_index: int) -> bool:
        return self._read_frame_info(physical_index).footer_validation in self._valid_footer_stamps

    def _frame_sort_key(self, physical_index: int) -> tuple[int, int, int]:
        frame = self._read_frame_info(physical_index)
        return frame.seconds, frame.subsecond_ticks, frame.first_record_number

    def _read_frame_info(self, physical_index: int) -> _FrameInfo:
        if not 0 <= physical_index < self.physical_frame_count:
            raise IndexError(f"Physical frame index {physical_index} is outside the file")
        frame_offset = (
            self.header.data_offset + physical_index * self.header.frame_size
        )
        with self.path.open("rb") as file:
            file.seek(frame_offset)
            header_bytes = file.read(_FRAME_HEADER_SIZE)
            file.seek(frame_offset + self.header.frame_size - _FRAME_FOOTER_SIZE)
            footer_bytes = file.read(_FRAME_FOOTER_SIZE)
        if len(header_bytes) != _FRAME_HEADER_SIZE or len(footer_bytes) != _FRAME_FOOTER_SIZE:
            raise ValueError(f"TOB3 frame {physical_index} is truncated")
        seconds, subsecond_ticks, first_record_number = _FRAME_HEADER.unpack(header_bytes)
        footer, footer_validation = _FRAME_FOOTER.unpack(footer_bytes)
        return _FrameInfo(
            seconds,
            subsecond_ticks,
            first_record_number,
            footer,
            footer_validation,
        )

    def _physical_frame_index(self, logical_frame_index: int) -> int:
        return (self._rotation + logical_frame_index) % self.frame_count

    def _timestamp_ns_at(self, record_index: int) -> int:
        if not 0 <= record_index < self.record_count:
            raise IndexError(f"Record index {record_index} is outside the file")
        logical_frame, in_frame = divmod(record_index, self.records_per_frame)
        frame = self._read_frame_info(self._physical_frame_index(logical_frame))
        return (
            int(_CAMPBELL_EPOCH_NS)
            + frame.seconds * 1_000_000_000
            + frame.subsecond_ticks * self._frame_timestamp_tick_ns
            + in_frame * self._record_interval_ns
        )

    def _time_lower_bound(self, value: object) -> int:
        target = _as_datetime64_ns(value).astype(np.int64)
        low, high = 0, self.record_count
        while low < high:
            middle = (low + high) // 2
            if self._timestamp_ns_at(middle) < target:
                low = middle + 1
            else:
                high = middle
        return low

    def _normalize_selection(
        self, *, records: slice | None, time: slice | None
    ) -> tuple[int, int]:
        if records is not None and time is not None:
            raise ValueError("Select by records or time, not both")
        if time is not None:
            if not isinstance(time, slice) or time.step is not None:
                raise ValueError("time must be a slice with no step")
            start = 0 if time.start is None else self._time_lower_bound(time.start)
            stop = self.record_count if time.stop is None else self._time_lower_bound(time.stop)
            return start, max(start, stop)

        if records is None:
            records = slice(None)
        if not isinstance(records, slice):
            raise TypeError("records must be a slice")
        start, stop, step = records.indices(self.record_count)
        if step != 1:
            raise ValueError("records slice step must be one")
        return start, stop

    def _select_fields(
        self, variables: Iterable[str] | None
    ) -> tuple[TOB3Field, ...]:
        if variables is None:
            return self.header.fields
        if isinstance(variables, str):
            variables = [variables]
        by_name = {field.name: field for field in self.header.fields}
        names = tuple(variables)
        unknown = [name for name in names if name not in by_name]
        if unknown:
            raise ValueError(f"Unknown TOB3 variables: {unknown}")
        return tuple(by_name[name] for name in names)

    def _frame_dtype(self) -> np.dtype:
        return np.dtype(
            {
                "names": [
                    "seconds",
                    "subsecond_ticks",
                    "first_record_number",
                    "data",
                    "footer",
                    "footer_validation",
                ],
                "formats": [
                    "<u4",
                    "<u4",
                    "<u4",
                    (self._frame_record_dtype, (self.records_per_frame,)),
                    "<u2",
                    "<u2",
                ],
                "offsets": [
                    0,
                    4,
                    8,
                    _FRAME_HEADER_SIZE,
                    self.header.frame_size - 4,
                    self.header.frame_size - 2,
                ],
                "itemsize": self.header.frame_size,
            }
        )

    def _read_physical_frames(
        self, physical_start: int, count: int, selected_fields: tuple[TOB3Field, ...]
    ) -> dict[str, np.ndarray]:
        byte_count = count * self.header.frame_size
        with self.path.open("rb") as file:
            file.seek(
                self.header.data_offset + physical_start * self.header.frame_size
            )
            raw = file.read(byte_count)
        if len(raw) != byte_count:
            raise ValueError("A selected TOB3 frame range is truncated")

        frames = np.frombuffer(raw, dtype=self._frame_dtype(), count=count)
        valid = np.isin(frames["footer_validation"], self._valid_footer_stamps)
        if not np.all(valid):
            bad = physical_start + int(np.flatnonzero(~valid)[0])
            raise ValueError(f"TOB3 frame {bad} has an invalid validation stamp")

        positions = np.arange(self.records_per_frame, dtype=np.uint64)
        record_numbers = (
            frames["first_record_number"].astype(np.uint64)[:, None] + positions
        ).astype(np.uint32, copy=False)
        timestamp_ns = (
            int(_CAMPBELL_EPOCH_NS)
            + frames["seconds"].astype(np.int64)[:, None] * 1_000_000_000
            + frames["subsecond_ticks"].astype(np.int64)[:, None]
            * self._frame_timestamp_tick_ns
            + positions.astype(np.int64)[None, :] * self._record_interval_ns
        )
        output = {
            "datetime": timestamp_ns.reshape(-1).view("datetime64[ns]").copy(),
            "record_number": record_numbers.reshape(-1).copy(),
        }
        for field in selected_fields:
            output[field.name] = _decode_field(
                frames["data"][field.name].reshape(-1), field.data_type
            )
        return output

    def _read_record_range(
        self, start: int, stop: int, selected_fields: tuple[TOB3Field, ...]
    ) -> dict[str, np.ndarray]:
        keys = ("datetime", "record_number") + tuple(
            field.name for field in selected_fields
        )
        if start == stop:
            return {
                "datetime": np.empty(0, dtype="datetime64[ns]"),
                "record_number": np.empty(0, dtype=np.uint32),
                **{
                    field.name: np.empty(0, dtype=_output_dtype(field.data_type))
                    for field in selected_fields
                },
            }

        first_logical_frame = start // self.records_per_frame
        last_logical_frame = (stop - 1) // self.records_per_frame
        frame_count = last_logical_frame - first_logical_frame + 1
        physical_start = self._physical_frame_index(first_logical_frame)
        first_run_count = min(frame_count, self.frame_count - physical_start)
        runs = [(physical_start, first_run_count)]
        if first_run_count < frame_count:
            runs.append((0, frame_count - first_run_count))

        run_arrays = [
            self._read_physical_frames(run_start, run_count, selected_fields)
            for run_start, run_count in runs
        ]
        combined = {
            key: (
                run_arrays[0][key]
                if len(run_arrays) == 1
                else np.concatenate([arrays[key] for arrays in run_arrays])
            )
            for key in keys
        }
        trim_start = start - first_logical_frame * self.records_per_frame
        trim_stop = trim_start + stop - start
        return {
            key: values[trim_start:trim_stop]
            for key, values in combined.items()
        }

    def _lazy_record_range(
        self,
        start: int,
        stop: int,
        selected_fields: tuple[TOB3Field, ...],
        chunks: int,
    ) -> dict[str, object]:
        # Keep the optional dask import out of header-only and eager reads.
        # atmPy's optional-import module also initializes its other scientific
        # facades, which can otherwise dominate the cost of inspecting a file.
        from atmPy import opt_imports

        dask = opt_imports.OptionalImport("dask", submodules="array")
        keys_and_dtypes = {
            "datetime": np.dtype("datetime64[ns]"),
            "record_number": np.dtype(np.uint32),
            **{
                field.name: _output_dtype(field.data_type)
                for field in selected_fields
            },
        }
        pieces = {key: [] for key in keys_and_dtypes}
        for chunk_start in range(start, stop, chunks):
            chunk_stop = min(chunk_start + chunks, stop)
            length = chunk_stop - chunk_start
            task = dask.delayed(self._read_record_range)(
                chunk_start, chunk_stop, selected_fields
            )
            for key, dtype in keys_and_dtypes.items():
                pieces[key].append(
                    dask.array.from_delayed(
                        dask.delayed(operator.getitem)(task, key),
                        shape=(length,),
                        dtype=dtype,
                    )
                )
        return {
            key: dask.array.concatenate(arrays)
            for key, arrays in pieces.items()
        }

    def _dataset_attrs(self) -> dict[str, object]:
        return {
            "source_file": str(self.path),
            "campbell_file_type": self.header.file_type,
            "station_name": self.header.station_name,
            "logger_model": self.header.logger_model,
            "logger_serial_number": self.header.logger_serial_number,
            "operating_system": self.header.operating_system,
            "program_name": self.header.program_name,
            "program_signature": self.header.program_signature,
            "file_created": self.header.file_created,
            "table_name": self.header.table_name,
            "record_interval": self.header.record_interval,
            "frame_time_resolution": self.header.frame_time_resolution,
            "frame_count": self.frame_count,
            "record_count": self.record_count,
            "record_order": "oldest_to_newest",
        }


def read_tob3(
    path: str | Path,
    *,
    records: slice | None = None,
    time: slice | None = None,
    variables: Iterable[str] | None = None,
    chunks: int | None = 250_000,
) -> xr.Dataset:
    """Read a Campbell TOB3 file or a chronological subset of it.

    This convenience function constructs :class:`TOB3Reader` and forwards all
    selection arguments to :meth:`TOB3Reader.read`.  See that method for the
    complete parameter and return-value documentation.
    """
    return TOB3Reader(path).read(
        records=records,
        time=time,
        variables=variables,
        chunks=chunks,
    )


def _record_dtype(fields: tuple[TOB3Field, ...]) -> np.dtype:
    return np.dtype(
        [(field.name, _storage_dtype(field.data_type)) for field in fields],
        align=False,
    )


def _dtype_with_itemsize(dtype: np.dtype, itemsize: int) -> np.dtype:
    """Return a structured dtype with trailing bytes added to each record."""
    if itemsize == dtype.itemsize:
        return dtype
    return np.dtype(
        {
            "names": list(dtype.names or ()),
            "formats": [dtype.fields[name][0] for name in dtype.names or ()],
            "offsets": [dtype.fields[name][1] for name in dtype.names or ()],
            "itemsize": itemsize,
        }
    )


def _storage_dtype(data_type: str) -> np.dtype:
    normalized = data_type.strip().upper()
    numeric = {
        "IEEE4": "<f4",
        "IEEE4L": "<f4",
        "IEEE4B": ">f4",
        "FP2": ">u2",
        "USHORT": "<u2",
        "SHORT": "<i2",
        "UINT2": ">u2",
        "INT2": ">i2",
        "ULONG": "<u4",
        "LONG": "<i4",
        "UINT4": ">u4",
        "INT4": ">i4",
        "BOOLEAN": "u1",
        "BOOL": "u1",
        "BOOL2": ">u2",
        "BOOL4": ">u4",
        "SECNANO": (">u4", (2,)),
        "NSEC": (">u4", (2,)),
    }
    if normalized in numeric:
        return np.dtype(numeric[normalized])
    match = _ASCII_TYPE.fullmatch(data_type.strip())
    if match:
        return np.dtype(f"S{int(match.group(1))}")
    raise NotImplementedError(f"Unsupported Campbell data type {data_type!r}")


def _output_dtype(data_type: str) -> np.dtype:
    normalized = data_type.strip().upper()
    if normalized == "FP2":
        return np.dtype(np.float32)
    if normalized in {"BOOLEAN", "BOOL", "BOOL2", "BOOL4"}:
        return np.dtype(bool)
    if normalized in {"SECNANO", "NSEC"}:
        return np.dtype("datetime64[ns]")
    match = _ASCII_TYPE.fullmatch(data_type.strip())
    if match:
        return np.dtype(f"U{int(match.group(1))}")
    return _storage_dtype(data_type).newbyteorder("=")


def _decode_field(values: np.ndarray, data_type: str) -> np.ndarray:
    normalized = data_type.strip().upper()
    if normalized == "FP2":
        return _decode_fp2(values)
    if normalized in {"BOOLEAN", "BOOL", "BOOL2", "BOOL4"}:
        return values.astype(bool, copy=True)
    if normalized in {"SECNANO", "NSEC"}:
        seconds = values[:, 0].astype(np.int64)
        nanoseconds = values[:, 1].astype(np.int64)
        timestamp_ns = int(_CAMPBELL_EPOCH_NS) + seconds * 1_000_000_000 + nanoseconds
        return timestamp_ns.view("datetime64[ns]")
    if _ASCII_TYPE.fullmatch(data_type.strip()):
        return np.char.rstrip(values.astype(_output_dtype(data_type)), "\x00")
    return values.astype(_output_dtype(data_type), copy=True)


def _decode_fp2(values: np.ndarray) -> np.ndarray:
    """Vectorized Campbell FP2 decoding, including INF and NAN sentinels."""
    raw = values.astype(np.uint16, copy=False)
    mantissa = raw & 0x1FFF
    exponent = (raw & 0x6000) >> 13
    decoded = mantissa.astype(np.float32) * np.power(
        np.float32(10.0), -exponent.astype(np.float32)
    )
    decoded = np.where(raw & 0x8000, -decoded, decoded).astype(np.float32)
    decoded[raw == 0x1FFF] = np.inf
    decoded[raw == 0x9FFF] = -np.inf
    decoded[raw == 0x9FFE] = np.nan
    return decoded


def _parse_record_interval(value: str) -> int:
    parts = value.strip().split()
    if len(parts) != 2:
        raise ValueError(f"Unsupported TOB3 record interval {value!r}")
    magnitude = decimal.Decimal(parts[0])
    unit_ns = {
        "NSEC": decimal.Decimal(1),
        "USEC": decimal.Decimal(1_000),
        "MSEC": decimal.Decimal(1_000_000),
        "SEC": decimal.Decimal(1_000_000_000),
        "MIN": decimal.Decimal(60_000_000_000),
        "HR": decimal.Decimal(3_600_000_000_000),
        "HOUR": decimal.Decimal(3_600_000_000_000),
        "DAY": decimal.Decimal(86_400_000_000_000),
    }
    try:
        interval_ns = magnitude * unit_ns[parts[1].upper()]
    except KeyError as error:
        raise ValueError(f"Unsupported TOB3 record interval unit {parts[1]!r}") from error
    if interval_ns != interval_ns.to_integral_value() or interval_ns <= 0:
        raise ValueError(f"TOB3 record interval {value!r} is not whole nanoseconds")
    return int(interval_ns)


def _parse_frame_resolution(value: str) -> int:
    match = _FRAME_RESOLUTION.fullmatch(value.strip())
    if not match:
        raise ValueError(f"Unsupported TOB3 frame time resolution {value!r}")
    scale = {"NSEC": 1, "USEC": 1_000, "MSEC": 1_000_000}
    magnitude = int(match.group(1) or 1)
    return magnitude * scale[match.group(2).upper()]


def _as_datetime64_ns(value: object) -> np.datetime64:
    if isinstance(value, dt.datetime) and value.tzinfo is not None:
        raise ValueError("TOB3 timestamps are naive; time slice bounds must not have tzinfo")
    try:
        return np.datetime64(value, "ns")
    except (TypeError, ValueError) as error:
        raise ValueError(f"Cannot interpret {value!r} as a datetime") from error


__all__ = ["TOB3Field", "TOB3Header", "TOB3Reader", "read_tob3", "read_tob3_header"]
