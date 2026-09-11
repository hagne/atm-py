"""
This module provides a class for managing a database of radflux parameters. It is likely that this class is not directly used by users,
but rather as a blueprint for a customized database class.
"""
import sqlite3
import socket
import pathlib as pl
import pandas as pd
import xarray as xr

radflux_parameter_table = 'radflux_parameters'
radflux_parameter_names = (
    'normalized_total_shortwave_power_coefficient',
    'normalized_total_shortwave_power_exponent',
    'normalized_diffuse_ratio_power_coefficient',
    'normalized_diffuse_ratio_power_exponent',
    'n_clear',
    'clear_fraction',
    'n_clear_global_irradiance_termporal_gradiant',
    'n_clear_normalized_diffuse_ratio_variability',
    'n_clear_diffuse_magnitude',
    'n_clear_normalized_global_magnitude',
    'mu0_coverage',
    'normalized_total_shortwave_median_absolute_deviation',
    'normalized_total_shortwave_coefficient_of_determination',
    'normalized_diffuse_ratio_median_absolute_deviation',
    'normalized_diffuse_ratio_coefficient_of_determination',
    'n_clear_above_limit',
    'diffuse_ratio_power_exponent_above_validity_limit',
    'diffuse_ratio_power_exponent_below_validity_limit',
    'normalized_total_shortwave_cosine_coverage_sufficient',
    'normalized_total_shortwave_final_iteration_cosine_coverage_sufficient',
    'diffuse_ratio_cosine_coverage_sufficient',
    'normalized_total_shortwave_power_coefficient_is_valid',
    'normalized_total_shortwave_power_exponent_is_valid',
    'normalized_diffuse_ratio_power_coefficient_is_valid',
    'normalized_diffuse_ratio_power_exponent_is_valid',
)
radflux_integer_parameter_names = frozenset({
    'n_clear',
    'n_clear_global_irradiance_termporal_gradiant',
    'n_clear_normalized_diffuse_ratio_variability',
    'n_clear_diffuse_magnitude',
    'n_clear_normalized_global_magnitude',
})
radflux_boolean_parameter_names = frozenset({
    'n_clear_above_limit',
    'diffuse_ratio_power_exponent_above_validity_limit',
    'diffuse_ratio_power_exponent_below_validity_limit',
    'normalized_total_shortwave_cosine_coverage_sufficient',
    'normalized_total_shortwave_final_iteration_cosine_coverage_sufficient',
    'diffuse_ratio_cosine_coverage_sufficient',
    'normalized_total_shortwave_power_coefficient_is_valid',
    'normalized_total_shortwave_power_exponent_is_valid',
    'normalized_diffuse_ratio_power_coefficient_is_valid',
    'normalized_diffuse_ratio_power_exponent_is_valid',
})
radflux_clearsky_parameter_names = frozenset({
    'normalized_diffuse_ratio_power_coefficient',
    'normalized_diffuse_ratio_power_exponent',
    'normalized_total_shortwave_power_coefficient',
    'normalized_total_shortwave_power_exponent'  
})
radflux_test_dict = {
    'normalized_total_shortwave_power_coefficient': 'normalized_total_shortwave_power_coefficient_is_valid',
    'normalized_total_shortwave_power_exponent':'normalized_total_shortwave_power_exponent_is_valid',
    'normalized_diffuse_ratio_power_coefficient': 'normalized_diffuse_ratio_power_coefficient_is_valid',
    'normalized_diffuse_ratio_power_exponent': 'normalized_diffuse_ratio_power_exponent_is_valid',}


radflux_float_parameter_names = frozenset(radflux_parameter_names).difference(
    radflux_integer_parameter_names,
    radflux_boolean_parameter_names,
)
radflux_discrete_parameter_names = (
    radflux_integer_parameter_names | radflux_boolean_parameter_names
)
radflux_table_columns = {
    'local_day': 'TEXT PRIMARY KEY',
    'input_file': 'TEXT NOT NULL',
    'next_day_needed': 'BOOLEAN',
    # 'output_file': 'TEXT NOT NULL',
    'processed_at': 'TEXT NOT NULL',
    'process_version': 'TEXT NOT NULL',
    'processing_server': 'TEXT NOT NULL',
    # 'clear_sky_params_optimized': 'TEXT',
    **{
        name: 'INTEGER' if name in radflux_discrete_parameter_names else 'REAL'
        for name in radflux_parameter_names
    },
}

class RadfluxParameterDatabase:
    def __init__(self, radflux_parameters_db, create_if_not_exist=False, version=None, verbose=False):
        self.radflux_parameters_db = pl.Path(radflux_parameters_db)
        self.verbose = verbose
        if version is None:
            version = '0.0'
        self.version = version
        if not self.radflux_parameters_db.exists() and not create_if_not_exist:
            raise FileNotFoundError(f"{self.radflux_parameters_db} does not exist and create_if_not_exist is False.")
        elif not self.radflux_parameters_db.exists() and create_if_not_exist:
            self.radflux_parameters_db.parent.mkdir(parents=True, exist_ok=True)
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)

    def connect_database(self):
        conn = sqlite3.connect(self.radflux_parameters_db)
        conn.row_factory = sqlite3.Row
        return conn

    def ensure_parameter_table(self, conn):
        parameter_columns = ',\n                '.join(
            f'{name} {dtype}' for name, dtype in radflux_table_columns.items()
        )
        table_exists = conn.execute(
            """
            SELECT name
            FROM sqlite_master
            WHERE type = 'table' AND name = ?
            """,
            (radflux_parameter_table,),
        ).fetchone() is not None
        conn.execute(f"""
            CREATE TABLE IF NOT EXISTS {radflux_parameter_table} (
                {parameter_columns}
            )
            """)
        existing_columns = {
            row['name']
            for row in conn.execute(f'PRAGMA table_info({radflux_parameter_table})')
        }
        if table_exists and 'local_day' not in existing_columns:
            if 'row_timestamp' not in existing_columns:
                raise ValueError(
                    f'{self.radflux_parameters_db} contains a '
                    f'{radflux_parameter_table} table without local_day'
                )
            conn.execute(
                f'ALTER TABLE {radflux_parameter_table} '
                'RENAME COLUMN row_timestamp TO local_day'
            )
            conn.execute(
                f'UPDATE {radflux_parameter_table} '
                'SET local_day = SUBSTR(local_day, 1, 10)'
            )
            existing_columns.remove('row_timestamp')
            existing_columns.add('local_day')
        for name, dtype in radflux_table_columns.items():
            if name not in existing_columns:
                dtype = dtype.replace(' NOT NULL', '').replace(' PRIMARY KEY', '')
                conn.execute(
                    f'ALTER TABLE {radflux_parameter_table} '
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
    
    def dump_radflux_parameters(self):
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)
            df = pd.read_sql_query(
                f'SELECT * FROM {radflux_parameter_table} ORDER BY local_day DESC',
                conn,
                index_col='local_day',
            )
        return df

    def get_clearsky_parameters(self, date, method = 'interpolate', ):
        

        local_day = self.date2dbformat(date)
        # selected_columns = ', '.join(('local_day', *radflux_parameter_names))
        out = {}
        for name in radflux_clearsky_parameter_names:
            selected_columns = f'local_day, {name}'
            test = radflux_test_dict[name]
            optimized_filter = (f"{test} IN ('True', 'true', 'TRUE', '1')")
            with self.connect_database() as conn:
                self.ensure_parameter_table(conn)
                previous = conn.execute(
                    f"""
                    SELECT {selected_columns}
                    FROM {radflux_parameter_table}
                    WHERE local_day <= ?
                    AND n_clear_above_limit IN ('True', 'true', 'TRUE', '1')
                    AND {optimized_filter}
                    ORDER BY local_day DESC
                    LIMIT 1
                    """,
                    (local_day,),
                ).fetchone()
                following = conn.execute(
                    f"""
                    SELECT {selected_columns}
                    FROM {radflux_parameter_table}
                    WHERE local_day >= ?
                    AND n_clear_above_limit IN ('True', 'true', 'TRUE', '1')
                    AND {optimized_filter}
                    ORDER BY local_day ASC
                    LIMIT 1
                    """,
                    (local_day,),
                ).fetchone()
            if method == 'previous':
                following = None
            elif method == 'following':
                previous = None
            elif method == 'interpolate':
                pass
            else:
                assert(False), f'Invalid method: {method}. Chose from "previous", "following", or "interpolate".'

            if previous is None and following is None:
                status = f'No optimized clearsky parameters found for {local_day}.'
                value = {name: None}       
            elif previous is None:
                value = following
                status = f'extrapolated, no previous parameters found, closest valid clearsky day: {following["local_day"]}'
            elif following is None:
                value = previous
                status = f'extrapolated, no following parameters found, closest valid clearsky day: {previous["local_day"]}'
            else:
                previous_time = pd.to_datetime(previous['local_day']).value
                following_time = pd.to_datetime(following['local_day']).value
                date_time = pd.to_datetime(date).value
                weight = (date_time - previous_time) / (following_time - previous_time)
                previous_value = previous[name]
                following_value = following[name]
                value = previous_value + (
                    following_value - previous_value
                ) * weight
                status = f'interpolated, closest valid clearsky days: {previous["local_day"]} and {following["local_day"]}'
            out[name] = {'value': value, 'status': status} 
        # return out
        ds = xr.Dataset()
        for k in out:
            ds[k] = xr.DataArray(out[k]['value'][k], attrs={'status': out[k]['status']})
        return ds


    # def read_previous_valid_clearsky_parameters(self, date):
    #     """Retrieves the last set of clearsky parameters before the given local day."""
    #     local_day = self.date2dbformat(date)
    #     optimized_filter = (
    #         # "clear_sky_params_optimized IN ('True', 'true', 'TRUE', '1')"
    #         "True"
    #     )
    #     with self.connect_database() as conn:
    #         self.ensure_parameter_table(conn)
    #         previous = conn.execute(
    #             f"""
    #             SELECT *
    #             FROM {radflux_parameter_table}
    #             WHERE local_day < ?
    #               AND {optimized_filter}
    #             ORDER BY local_day DESC
    #             LIMIT 1
    #             """,
    #             (local_day,),
    #         ).fetchone()
    #         self.tp_prvious = previous

    #     if previous is None:
    #         if self.verbose:
    #             print(f'No previous optimized clearsky parameters found for {local_day}.')
    #         self.tp_previous_radflux_parameters_record = None
    #         return None

    #     self.tp_previous_radflux_parameters_record = dict(previous)
    #     return self._optimization_results_dataset(
    #         previous,
    #         f'previous valid clearsky day: {previous["local_day"]}',
    #     )
    
    def write_radflux_parameters(self,
                                  date,
                                  path2file,
                                  clearsky_parameters,
                                  processing_date=None,
                                  processing_server=None,
                                  next_day_needed=None
                                  ):

        """Insert or update one set of radflux processing parameters.

        Parameters
        ----------
        date : datetime-like
            Local day covered by the processed data.
        path2file : path-like
            Path to the processed input file.
        clearsky_parameters : xarray.Dataset or None
            Clear-sky optimization results. The dataset must contain every
            variable named in ``radflux_parameter_names``. If ``None``, the
            corresponding database columns are stored as null values and the
            day is marked as not optimized.
        processing_date : str, optional
            Timestamp describing when the record was processed. Defaults to
            the current timestamp.
        processing_server : str, optional
            Name of the server that processed the record. Defaults to the
            current host name.
        next_day_needed : bool
            Whether processing requires data from the following day.

        Notes
        -----
        An existing record for the same local day is replaced in place.
        """

        if processing_date is None:
            processing_date = pd.Timestamp.now().isoformat()
        if processing_server is None:
            processing_server = socket.gethostname()

        if clearsky_parameters is None:
            parameters = {}
        elif not isinstance(clearsky_parameters, xr.Dataset):
            raise TypeError('clearsky_parameters must be an xarray.Dataset or None')
        else:
            missing = set(radflux_parameter_names).difference(
                clearsky_parameters.data_vars
            )
            if missing:
                raise ValueError(
                    'clearsky_parameters is missing optimization results: '
                    f'{sorted(missing)}'
                )
            parameters = {
                name: self._database_value(clearsky_parameters[name])
                for name in radflux_parameter_names
            }
        values = {
            'local_day': self.date2dbformat(date),
            'input_file': str(path2file),
            'next_day_needed': next_day_needed,
            # 'output_file': str(row.p2f_out),
            'processed_at': processing_date,
            'process_version': self.version,
            'processing_server': processing_server,
            # 'clear_sky_params_optimized': clearsky_parameters is not None,
            # 'parameters_json': json.dumps(parameters, sort_keys=True),
        }
        values.update({
            name: parameters.get(name)
            for name in radflux_parameter_names
        })
        columns = tuple(values)
        placeholders = ', '.join(['?'] * len(columns))
        update_columns = ', '.join(
            f'{column} = excluded.{column}'
            for column in columns
            if column != 'local_day'
        )

        self.tp_columns = columns
        self.tp_placeholders = placeholders
        self.tp_values = values
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)
            conn.execute(
                f"""
                INSERT INTO {radflux_parameter_table}
                    ({', '.join(columns)})
                VALUES ({placeholders})
                ON CONFLICT(local_day) DO UPDATE SET
                    {update_columns}
                """,
                tuple(values[column] for column in columns),
            )
    def delete_rows_on_date(self, date, find_only = False):
        """Deletes all rows in the radflux parameter database for a specific date."""
        date_str = pd.to_datetime(date).strftime('%Y-%m-%d')
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)
            if find_only:
                rows = conn.execute(
                    f"""
                    SELECT *
                    FROM {radflux_parameter_table}
                    WHERE local_day = ?
                    """,
                    (date_str,),
                ).fetchall()
                return [dict(row) for row in rows]
            else:
                conn.execute(
                    f"""
                    DELETE FROM {radflux_parameter_table}
                    WHERE local_day = ?
                    """,
                    (date_str,),
                )
        
