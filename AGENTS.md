# AGENTS.md

Guidance for future edits in `atmPy`.

## Local Style First

Before editing, inspect nearby code and follow the dominant local pattern.

Prefer direct, idiomatic xarray expressions over generic helper functions.

Prefer:

```python
aod = ds.aod.sel(datetime=time, wavelength=wl)
```

Avoid:

```python
aod = get_aod(ds, time=time, wavelength=wl)
```

unless the helper already exists or the operation is reused and nontrivial.

Do not add new abstractions, wrappers, compatibility layers, or broad refactors unless explicitly requested.

## Edit Scope

- Keep changes local to the requested module or behavior.
- Do not reformat old files, remove old comments, or migrate broad style unless asked.
- When touching old code, flag nearby outdated patterns in the final response instead of silently rewriting them.

Example flag:

```text
Note: this function still converts an xarray Dataset to pandas before grouping.
Future work should keep this in xarray so dask chunks remain lazy.
```

## Python And Imports

- Target the Python baseline in `pyproject.toml` (`>=3.11`). New code may use modern typing such as `list[str]`, `str | None`, and `typing.Self`.
- Match the import style already used in the file being edited. Do not introduce repo-wide import alias changes unless explicitly asked.
- Prefer keeping imported submodules in their parent package namespace.

Prefer:

```python
import matplotlib as mpl
import matplotlib.colors

norm = mpl.colors.Normalize(vmin=0, vmax=1)
```

Avoid:

```python
from matplotlib import colors

norm = colors.Normalize(vmin=0, vmax=1)
```

- Prefer `pathlib.Path` for new path handling. Use `os` when it is the right interface, for example `os.environ`.
- Preserve existing public APIs. If replacing legacy behavior, keep a clear compatibility wrapper or explicit deprecation error.

## Data Model

- Store, operate on, and return scientific data as `xarray.Dataset` or `xarray.DataArray`.
- Prefer direct xarray access and selection over helper accessors.
- Preserve dims, coords, attrs, units, and dask-backed laziness.
- Do not add flexible variable-name helpers or functions that accept alternate names for the same scientific variable. Use the exact expected variable names in the dataset.
- Validate units attributes explicitly. If a variable exists but its `units` attribute does not match the expected units, raise an error rather than converting between units.
- Use `pandas` only when `xarray` lacks the operation or when a measured performance issue justifies the conversion.
- Use raw `numpy` arrays only for low-level kernels or library APIs that require `ndarray`; wrap kernels with `xarray.apply_ufunc` or `xarray.map_blocks` when possible.

Prefer:

```python
ds = _xr.Dataset(
    data_vars={"aod": (("datetime", "wavelength"), values)},
    coords={"datetime": times, "wavelength": wavelengths},
    attrs={"source": "surfrad"},
)
daily = ds.resample(datetime="1D").mean()
```

Avoid:

```python
df = ds.to_dataframe()
daily = df.groupby(df.index.date).mean()
return daily.values
```

unless the conversion is necessary and documented.

## Large Data

- Prefer chunk-aware IO: `xarray.open_dataset(..., chunks=...)`, `xarray.open_mfdataset(...)`, NetCDF, and Zarr.
- Avoid accidental eager loads in core operations: `.values`, `.to_numpy()`, `.compute()`, `.load()`, `_np.asarray(...)`, `list(...)`, or loops over large dimensions.

Prefer:

```python
ds = _xr.open_mfdataset(paths, chunks={"datetime": "auto"})
out = ds.rolling(datetime=24).mean()
```

For custom kernels:

```python
out = _xr.apply_ufunc(
    kernel,
    ds["aod"],
    input_core_dims=[["wavelength"]],
    output_core_dims=[["wavelength"]],
    dask="parallelized",
)
```

## Optional Dependencies

- Follow the project’s existing optional-dependency pattern.
- If the project has an `opt_imports` module, use it for optional scientific packages so imports fail only when the feature is used.
- When optional dependencies need submodules, load them through the optional-import mechanism rather than importing submodules directly.
- For new optional feature groups, add a matching `pyproject.toml` optional dependency group when appropriate.

Prefer predefined optional imports when they already exist:

```python
from atmPy.opt_imports import plt

plt.plot(x, y)
```

Prefer explicit OptionalImport objects for optional packages or submodules that do not yet have predefined imports:

```python
from atmPy import opt_imports

geopy = opt_imports.OptionalImport("geopy", submodules="distance")

d = geopy.distance.distance(point_a, point_b)
```

Avoid direct submodule imports for optional dependencies:

```python
from geopy import distance

d = distance.distance(point_a, point_b)
```

## Cached Data

- Keep expensive loads/calculations cached behind properties.
- Invalidate cached data in setters that change inputs.
- Cache private intermediate objects on private attributes.
- Cache derived scientific variables in the `xarray.Dataset` when they belong with the dataset.

Prefer private-attribute caching for loaded objects:

```python
@property
def aod(self):
    if self._aod is None:
        self._aod = self.load_data(param="aod")
    return self._aod
```

Prefer dataset caching for derived scientific variables:

```python
@property
def angstrom_exponent(self):
    if "angstrom_exponent" not in self.ds:
        self.ds["angstrom_exponent"] = self._calc_angstrom_exponent()
    return self.ds.angstrom_exponent
```

## Tests

- Add or update focused `unittest` tests under `tests/` for changed calculations.
- Treat expected numerical values as regression baselines. Do not update expected numbers just to make tests pass after a code change.
- Change established expected values only when the scientific/algorithmic change is deliberate, documented, and approved by the maintainer.
- Prioritize tests that verify scientific results and numerical data values over tests that only verify structure, metadata, or implementation details.
- Do not add many trivial asserts. Prefer a few meaningful assertions that test the scientific result, numerical behavior, and key coordinates/metadata only when relevant.
- Run:

```bash
./run_tests.sh
```

- Stub optional dependencies in tests when a unit test does not need the real package.

## Final Response

Keep final responses concise:
- summarize what changed
- mention tests run
- flag nearby outdated patterns only when relevant