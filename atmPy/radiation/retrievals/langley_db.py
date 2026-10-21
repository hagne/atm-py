"""SQLite storage for Langley calibration results."""

import pathlib as pl
import socket
import sqlite3
from contextlib import closing

import pandas as pd
import xarray as xr


langley_result_table = 'langley_results'
langley_fitresult_names = (
    'slope',
    'intercept',
    'slope_stderr',
    'intercept_stderr',
    'resid_curvature',
)
airmass_coverage_metric_names = (
    'retained_points',
    'retained_fraction',
    'cdf_span',
    'maximum_cdf_gap',
    'maximum_cdf_gap_ratio',
    'cdf_spacing_cv',
    'effective_airmass',
    'effective_airmass_offset',
)
langley_integer_result_names = frozenset({
    'retained_points',
})
langley_result_names = (
    *langley_fitresult_names,
    *airmass_coverage_metric_names,
    'residual_correlation_determinant',
)
langley_table_columns = {
    'local_day': 'TEXT NOT NULL',
    'time_of_day': 'TEXT NOT NULL',
    'channel': 'REAL NOT NULL',
    'input_file': 'TEXT NOT NULL',
    'processed_at': 'TEXT NOT NULL',
    'process_version': 'TEXT NOT NULL',
    'processing_server': 'TEXT NOT NULL',
    **{
        name: 'INTEGER' if name in langley_integer_result_names else 'REAL'
        for name in langley_result_names
    },
}

class LangleyResultDatabase:
    """Manage channel-resolved Langley results in a SQLite database.

    Each database row contains the fit results and airmass-coverage metrics
    for one local day, time of day, and instrument channel. Writing the same
    combination of ``local_day``, ``time_of_day``, and ``channel`` again
    updates the existing row.

    Parameters
    ----------
    langley_results_db : path-like
        Path to the SQLite database file.
    create_if_not_exist : bool, optional
        Create the database and its parent directory when they do not exist.
        The default is ``False``.
    version : str, optional
        Processing-version identifier stored with every result. The default
        is ``"0.0"``.
    verbose : bool, optional
        Enable informational output. Reserved for compatibility with the
        other database interfaces; the default is ``False``.

    Examples
    --------
    Given an existing :class:`Langley` instance named ``langley``, create a
    database, store its results, and inspect the stored rows:

    >>> database = LangleyResultDatabase(
    ...     'langley_results.sqlite',
    ...     create_if_not_exist=True,
    ...     version='1.0',
    ... )
    >>> database.write_langley(
    ...     date='2026-05-02',
    ...     path2file='mfrsr_20260502.nc',
    ...     langley=langley,
    ... )
    >>> results = database.dump_langley_results()
    >>> results.loc[('2026-05-02', 'am')]
    """

    def __init__(self, langley_results_db, create_if_not_exist=False,
                 version=None, verbose=False):
        self.langley_results_db = pl.Path(langley_results_db)
        self.verbose = verbose
        if version is None:
            version = '0.0'
        self.version = version
        if not self.langley_results_db.exists() and not create_if_not_exist:
            raise FileNotFoundError(
                f'{self.langley_results_db} does not exist and '
                'create_if_not_exist is False.'
            )
        if not self.langley_results_db.exists():
            self.langley_results_db.parent.mkdir(parents=True, exist_ok=True)
        with closing(self.connect_database()) as conn:
            with conn:
                self.ensure_result_table(conn)

    def connect_database(self):
        conn = sqlite3.connect(self.langley_results_db)
        conn.row_factory = sqlite3.Row
        return conn

    def ensure_result_table(self, conn):
        result_columns = ',\n                '.join(
            f'{name} {dtype}' for name, dtype in langley_table_columns.items()
        )
        conn.execute(f"""
            CREATE TABLE IF NOT EXISTS {langley_result_table} (
                {result_columns},
                PRIMARY KEY (local_day, time_of_day, channel)
            )
            """)
        existing_columns = {
            row['name']
            for row in conn.execute(f'PRAGMA table_info({langley_result_table})')
        }
        for name, dtype in langley_table_columns.items():
            if name not in existing_columns:
                dtype = dtype.replace(' NOT NULL', '')
                conn.execute(
                    f'ALTER TABLE {langley_result_table} '
                    f'ADD COLUMN {name} {dtype}'
                )

    @staticmethod
    def date2dbformat(date):
        return pd.to_datetime(date).date().isoformat()

    @staticmethod
    def _database_value(value):
        if hasattr(value, 'item'):
            value = value.item()
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            pass
        return value

    @staticmethod
    def _validate_results(langley_fitresults, airmass_coverage_metrics):
        if not isinstance(langley_fitresults, pd.DataFrame):
            raise TypeError('langley_fitresults must be a pandas.DataFrame')
        missing_fitresults = set(langley_fitresult_names).difference(
            langley_fitresults.columns
        )
        if missing_fitresults:
            raise ValueError(
                'langley_fitresults is missing results: '
                f'{sorted(missing_fitresults)}'
            )
        if langley_fitresults.index.has_duplicates:
            raise ValueError('langley_fitresults contains duplicate channels')

        if not isinstance(airmass_coverage_metrics, xr.Dataset):
            raise TypeError('airmass_coverage_metrics must be an xarray.Dataset')
        missing_metrics = set(airmass_coverage_metric_names).difference(
            airmass_coverage_metrics.data_vars
        )
        if missing_metrics:
            raise ValueError(
                'airmass_coverage_metrics is missing metrics: '
                f'{sorted(missing_metrics)}'
            )
        if 'channel' not in airmass_coverage_metrics.coords:
            raise ValueError('airmass_coverage_metrics must have a channel coordinate')

        fit_channels = pd.Index(langley_fitresults.index)
        metric_channels = airmass_coverage_metrics.channel.to_index()
        if not fit_channels.equals(metric_channels):
            raise ValueError(
                'langley_fitresults and airmass_coverage_metrics must have '
                'the same channels in the same order'
            )

    def dump_langley_results(self, local_day_range=None, channel=None):
        """Return stored results, optionally filtered, newest day first.

        Parameters
        ----------
        local_day_range : tuple, optional
            Inclusive start and end dates.
        channel : float, optional
            A single channel to return.
        """
        where = []
        parameters = []
        if local_day_range is not None:
            start, end = local_day_range
            where.append('local_day BETWEEN ? AND ?')
            parameters.extend([
                self.date2dbformat(start),
                self.date2dbformat(end),
            ])
        if channel is not None:
            where.append('channel = ?')
            parameters.append(self._database_value(channel))
        where_clause = f"WHERE {' AND '.join(where)}" if where else ''

        with closing(self.connect_database()) as conn:
            with conn:
                self.ensure_result_table(conn)
                results = pd.read_sql_query(
                    f"""
                    SELECT *
                    FROM {langley_result_table}
                    {where_clause}
                    ORDER BY local_day DESC, time_of_day, channel
                    """,
                    conn,
                    params=parameters,
                )
        return results#.set_index(['local_day', 'time_of_day', 'channel'])

    def write_langley(self, date, path2file, langley, processing_date=None,
                      processing_server=None):
        """Insert or update the fit and coverage results from a Langley object."""
        if langley.when is None:
            raise ValueError('langley.when must identify the time of day')
        return self.write_langley_results(
            date=date,
            path2file=path2file,
            when=langley.when,
            langley_fitresults=langley.langley_fitres,
            airmass_coverage_metrics=langley.airmass_coverage_metrics,
            residual_correlation_determinant=(
                langley.langley_residual_correlation_prop['determinant']
            ),
            processing_date=processing_date,
            processing_server=processing_server,
        )

    def add_empty_entry(self, date, path2file, when, channel,
                        processing_date=None, processing_server=None):
        """Store a missing-result marker for one input file and channel."""
        if not isinstance(when, str) or not when:
            raise ValueError('when must be a non-empty string')
        if processing_date is None:
            processing_date = pd.Timestamp.now().isoformat()
        if processing_server is None:
            processing_server = socket.gethostname()

        key_columns = ('local_day', 'time_of_day', 'channel')
        provenance_columns = (
            'input_file',
            'processed_at',
            'process_version',
            'processing_server',
        )
        update_values = ', '.join((
            *(
                f'{column} = excluded.{column}'
                for column in provenance_columns
            ),
            *(f'{column} = NULL' for column in langley_result_names),
        ))
        with closing(self.connect_database()) as conn:
            with conn:
                self.ensure_result_table(conn)
                conn.execute(
                    f"""
                    INSERT INTO {langley_result_table}
                        ({', '.join((*key_columns, *provenance_columns))})
                    VALUES (?, ?, ?, ?, ?, ?, ?)
                    ON CONFLICT(local_day, time_of_day, channel) DO UPDATE SET
                        {update_values}
                    """,
                    (
                        self.date2dbformat(date),
                        when,
                        self._database_value(channel),
                        str(path2file),
                        processing_date,
                        self.version,
                        processing_server,
                    ),
                )

    def write_langley_results(self, date, path2file, when,
                              langley_fitresults, airmass_coverage_metrics,
                              residual_correlation_determinant=None,
                              processing_date=None, processing_server=None):
        """Insert or update one AM or PM set of channel-resolved results.

        Parameters
        ----------
        date : datetime-like
            Local day covered by the Langley analysis.
        path2file : path-like
            Path to the processed input file.
        when : str
            Time-of-day identifier, normally ``"am"`` or ``"pm"``.
        langley_fitresults : pandas.DataFrame
            Fit results indexed by channel, as returned by
            :attr:`Langley.langley_fitres`.
        airmass_coverage_metrics : xarray.Dataset
            Channel-resolved metrics returned by
            :attr:`Langley.airmass_coverage_metrics`.
        residual_correlation_determinant : float, optional
            Determinant of the cross-channel residual correlation matrix.
        processing_date : str, optional
            Processing timestamp. Defaults to the current timestamp.
        processing_server : str, optional
            Processing host. Defaults to the current host name.

        Notes
        -----
        Rows are updated in place for matching local day, time of day, and
        channel.
        """
        if not isinstance(when, str) or not when:
            raise ValueError('when must be a non-empty string')
        self._validate_results(langley_fitresults, airmass_coverage_metrics)
        if processing_date is None:
            processing_date = pd.Timestamp.now().isoformat()
        if processing_server is None:
            processing_server = socket.gethostname()

        common_values = {
            'local_day': self.date2dbformat(date),
            'time_of_day': when,
            'input_file': str(path2file),
            'processed_at': processing_date,
            'process_version': self.version,
            'processing_server': processing_server,
            'residual_correlation_determinant': self._database_value(
                residual_correlation_determinant
            ),
        }
        rows = []
        for channel, fitresults in langley_fitresults.iterrows():
            coverage = airmass_coverage_metrics.sel(channel=channel)
            values = dict(common_values)
            values['channel'] = self._database_value(channel)
            values.update({
                name: self._database_value(fitresults[name])
                for name in langley_fitresult_names
            })
            values.update({
                name: self._database_value(coverage[name])
                for name in airmass_coverage_metric_names
            })
            rows.append(values)

        columns = tuple(rows[0]) if rows else tuple(langley_table_columns)
        placeholders = ', '.join(['?'] * len(columns))
        update_columns = ', '.join(
            f'{column} = excluded.{column}'
            for column in columns
            if column not in {'local_day', 'time_of_day', 'channel'}
        )
        with closing(self.connect_database()) as conn:
            with conn:
                self.ensure_result_table(conn)
                conn.executemany(
                    f"""
                    INSERT INTO {langley_result_table}
                        ({', '.join(columns)})
                    VALUES ({placeholders})
                    ON CONFLICT(local_day, time_of_day, channel) DO UPDATE SET
                        {update_columns}
                    """,
                    [tuple(row[column] for column in columns) for row in rows],
                )

    def delete_rows_on_date(self, date, when=None, find_only=False):
        """Find or delete Langley rows for one local day and optional period."""
        local_day = self.date2dbformat(date)
        where = 'local_day = ?'
        parameters = [local_day]
        if when is not None:
            where += ' AND time_of_day = ?'
            parameters.append(when)

        with closing(self.connect_database()) as conn:
            with conn:
                self.ensure_result_table(conn)
                if find_only:
                    rows = conn.execute(
                        f'SELECT * FROM {langley_result_table} WHERE {where}',
                        parameters,
                    ).fetchall()
                    return [dict(row) for row in rows]
                conn.execute(
                    f'DELETE FROM {langley_result_table} WHERE {where}',
                    parameters,
                )
