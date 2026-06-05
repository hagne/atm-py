Radiation
---------

The :mod:`atmPy.radiation` package contains solar-position helpers,
shortwave-radiation retrievals, radiometer instrumentation helpers, Rayleigh
and Mie scattering calculations, and surface-albedo utilities.

Tilt Correction
===============

Moving platforms measure global shortwave irradiance on a tilted detector.
The tilt correction implemented in
:mod:`atmPy.radiation.retrievals.tiltcorrection` follows Long et al. (2010) for
downwelling shortwave irradiance. It assumes the direct beam is corrected by
the projection factor of the tilted detector while diffuse irradiance is
isotropic and therefore unaffected by tilt.

For a tilted global measurement :math:`G_t`, diffuse fraction :math:`f_D`,
horizontal projection factor :math:`\mu_0`, and tilted projection factor
:math:`\mu_t`, the corrected terms are:

.. math::

   G_0 = G_t \frac{-f_D\mu_0 + f_D\mu_t + \mu_0}{\mu_t}

.. math::

   N = G_t \frac{1 - f_D}{\mu_t}

.. math::

   D = G_t f_D

where :math:`G_0` is corrected global horizontal irradiance, :math:`N` is
direct normal irradiance, and :math:`D` is diffuse horizontal irradiance.

The input dataset must contain the exact variables used by the retrieval -- NOTE: solar position variables are provided by radiation objects when tilt correction is applied in the wrapper approach (first code block below):

.. list-table::
   :header-rows: 1

   * - Variable
     - Units
   * - ``global_horizontal``
     - ``W/m^2``
   * - ``diffuse_horizontal``
     - ``W/m^2``
   * - ``solar_zenith``
     - ``radian``
   * - ``solar_azimuth``
     - ``radian``
   * - ``platform_pitch``
     - ``radian``
   * - ``platform_roll``
     - ``radian``
   * - ``platform_heading``
     - ``radian``

Use the broadband-radiation wrapper when working with a combined global,
diffuse, and direct dataset:

.. code-block:: python

   from atmPy.radiation.retrievals import broadband_shortwave_radiation as bsr

   radiation = bsr.CombinedGlobalDiffuseDirect(ds)
   corrected = radiation.apply_tilt_correction(sensor_response_time=rs)

Or call the retrieval function directly:

.. code-block:: python

   from atmPy.radiation.retrievals import tiltcorrection

   corrected = tiltcorrection.apply_tilt_correction(
       ds,
       sensor_response_time=rs,
   )

The complete derivation and usage notes are also available as a notebook:

.. toctree::
   :maxdepth: 1

   notebooks/radiation/tilt_correction.ipynb

Module Overview
===============

.. toctree::
   :maxdepth: 2

   radiation_api.rst

.. automodule:: atmPy.radiation
   :members:
   :undoc-members:
   :show-inheritance:
