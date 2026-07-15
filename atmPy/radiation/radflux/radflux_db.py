"""
This module provides a class for managing a database of radflux parameters. It is likely that this class is not directly used by users,
but rather as a blueprint for a customized database class.
"""
import sqlite3
import pathlib as pl
import pandas as pd


default_clearsky_params = {'nsw_exp': 1.202095545434091,
                            'nsw_min': 800,
                            'nsw_max': 1400,
                            'ndr_exp': -0.6827046137686424,
                            'mu0_min': 0.05,
                            'diffuse_max_coeff': 150,
                            'diffuse_max_exp': 0.5,
                            'max_dsw_dt': 8,
                            'ndr_std_max': 0.005,
                            'ndr_window': 11,}


radflux_parameter_table = 'radflux_parameters'
radflux_parameter_names = tuple(default_clearsky_params) + (
    'nsw_coeff',
    'nsw_r2',
    'ndr_std_max_estimated',
    # 'diffuse_max_coeff_estimated',
    # 'diffuse_max_exp_estimated',
    'normalized_diffuse_fit_exp',
    'normalized_diffuse_fit_coeff',
    'max_dsw_dt_estimated',
)
radflux_table_columns = {
    'row_timestamp': 'TEXT PRIMARY KEY',
    'input_file': 'TEXT NOT NULL',
    'next_day_needed': 'BOOLEAN',
    # 'output_file': 'TEXT NOT NULL',
    'processed_at': 'TEXT NOT NULL',
    'process_version': 'TEXT NOT NULL',
    'processing_server': 'TEXT NOT NULL',
    'clear_sky_params_optimized': 'TEXT',
    # 'parameters_json': 'TEXT NOT NULL',
    **{name: 'REAL' for name in radflux_parameter_names},
}

class RadfluxParameterDatabase:
    def __init__(self, radflux_parameters_db, create_if_not_exist=False, version=None, verbose=False):
        self.radflux_parameters_db = pl.Path(radflux_parameters_db)
        self.verbose = verbose
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
        if table_exists and 'row_timestamp' not in existing_columns:
            raise ValueError(
                f'{self.radflux_parameters_db} contains a '
                f'{radflux_parameter_table} table without row_timestamp'
            )
        for name, dtype in radflux_table_columns.items():
            if name not in existing_columns:
                dtype = dtype.replace(' NOT NULL', '').replace(' PRIMARY KEY', '')
                conn.execute(
                    f'ALTER TABLE {radflux_parameter_table} '
                    f'ADD COLUMN {name} {dtype}'
                )

    @staticmethod
    def timestamp2dbformat(timestamp):
        return pd.to_datetime(timestamp).isoformat()

    @staticmethod
    def _database_value(value):
        try:
            if pd.isna(value):
                return None
        except (TypeError, ValueError):
            pass
        if hasattr(value, 'item'):
            value = value.item()
        return value
    
    def dump_radflux_parameters(self):
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)
            df = pd.read_sql_query(
                f'SELECT * FROM {radflux_parameter_table} ORDER BY row_timestamp DESC',
                conn,
                index_col='row_timestamp',
            )
        return df

    def get_clearsky_parameter(self, date):
        """Returns clearsky parameters for the given day. Parmeters are either interpolated, 
        when a valid paremeters exist before and after the date, or extrapolated (same as 
        last/first valid paremeters), when only one (before/after) valid set of parameters exists.
        
        Parameters
        ----------
        date: datetime-like
            The date for which to retrieve the clearsky parameters.

        Returns
        -------
        dict
            A dictionary containing the clearsky parameters for the specified date. This includes 
            an entry 'status' which indicates whether the parameters were interpolated, extrapolated, 
            or if no valid parameters were found.
        """
        row_timestamp = self.timestamp2dbformat(date)
        selected_columns = ', '.join(('row_timestamp', *radflux_parameter_names))
        optimized_filter = (
            "clear_sky_params_optimized IN ('True', 'true', 'TRUE', '1')"
        )
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)
            previous = conn.execute(
                f"""
                SELECT {selected_columns}
                FROM {radflux_parameter_table}
                WHERE row_timestamp <= ?
                  AND {optimized_filter}
                ORDER BY row_timestamp DESC
                LIMIT 1
                """,
                (row_timestamp,),
            ).fetchone()
            following = conn.execute(
                f"""
                SELECT {selected_columns}
                FROM {radflux_parameter_table}
                WHERE row_timestamp >= ?
                  AND {optimized_filter}
                ORDER BY row_timestamp ASC
                LIMIT 1
                """,
                (row_timestamp,),
            ).fetchone()

        if previous is None and following is None:
            if self.verbose:
                print(f'No optimized clearsky parameters found for {row_timestamp}.')
            parameters = default_clearsky_params.copy()
            parameters['status'] = 'no valid parameters found'
            return parameters

        def row_to_parameters(row, status):
            # parameters = default_clearsky_params.copy()
            parameters = {
                name: row[name]
                for name in radflux_parameter_names
                if row[name] is not None
            }
            parameters['status'] = status
            return parameters

        if previous is None:
            return row_to_parameters(following, f'extrapolated, no previous parameters found, closest valid clearsky day: {following["row_timestamp"]}')
        if following is None:
            return row_to_parameters(previous, f'extrapolated, no following parameters found, closest valid clearsky day: {previous["row_timestamp"]}')
        if previous['row_timestamp'] == following['row_timestamp']:
            return row_to_parameters(previous, 'valid clearsky day, no interpolation needed')

        previous_time = pd.to_datetime(previous['row_timestamp']).value
        following_time = pd.to_datetime(following['row_timestamp']).value
        date_time = pd.to_datetime(date).value
        weight = (date_time - previous_time) / (following_time - previous_time)

        parameters = default_clearsky_params.copy()
        for name in radflux_parameter_names:
            previous_value = previous[name]
            following_value = following[name]
            if previous_value is None and following_value is None:
                continue
            if previous_value is None:
                parameters[name] = following_value
            elif following_value is None:
                parameters[name] = previous_value
            else:
                parameters[name] = previous_value + (
                    following_value - previous_value
                ) * weight
        parameters['status'] = f'interpolated, closest valid clearsky days: {previous["row_timestamp"]} and {following["row_timestamp"]}'
        return parameters

    def read_previous_valid_clearsky_parameters(self, timestamp):
        """Retrieves the last set of clearsky parameters before the given timestamp."""
        row_timestamp = self.timestamp2dbformat(timestamp)
        with self.connect_database() as conn:
            self.ensure_parameter_table(conn)
            previous = conn.execute(
                f"""
                SELECT *
                FROM {radflux_parameter_table}
                WHERE row_timestamp < ?
                  AND clear_sky_params_optimized = 'True'
                ORDER BY row_timestamp DESC
                LIMIT 1
                """,
                (row_timestamp,),
            ).fetchone()
            self.tp_prvious = previous

        if previous is None:
            if self.verbose:
                print(f'No previous optimized clearsky parameters found for {row_timestamp}.')
            self.tp_previous_radflux_parameters_record = None
            return None

        self.tp_previous_radflux_parameters_record = dict(previous)
        # parameters = json.loads(previous['parameters_json'])
        # parameters = {
        #     name: value
        #     for name, value in parameters.items()
        #     if value is not None
        # }
        # return {**default_clearsky_params, **parameters}
        return dict(previous)

    def write_radflux_parameters(self,
                                  row,
                                  clearsky_parameters,
                                  processing_date,
                                  processing_server,
                                  clear_sky_params_optimized,
                                  next_day_needed
                                  ):
        # cleanup the parameters to ensure they are JSON serializable and handle NaN values
        parameters = {
            name: self._database_value(value)
            for name, value in clearsky_parameters.items()
        }
        values = {
            'row_timestamp': self.timestamp2dbformat(row.name),
            'input_file': str(row.p2f_in),
            'next_day_needed': next_day_needed,
            # 'output_file': str(row.p2f_out),
            'processed_at': processing_date,
            'process_version': self.version,
            'processing_server': processing_server,
            'clear_sky_params_optimized': clear_sky_params_optimized,
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
            if column != 'row_timestamp'
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
                ON CONFLICT(row_timestamp) DO UPDATE SET
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
                    WHERE DATE(row_timestamp) = ?
                    """,
                    (date_str,),
                ).fetchall()
                return [dict(row) for row in rows]
            else:
                conn.execute(
                    f"""
                    DELETE FROM {radflux_parameter_table}
                    WHERE DATE(row_timestamp) = ?
                    """,
                    (date_str,),
                )
        
