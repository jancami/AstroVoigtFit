master_function
===============

.. autofunction:: astrovoigtfit.master_function

This function is a wrapper around ``mother_function`` and is optimized for distributing
absorption line parameters across multiple species and components (clouds) in a flexible,
vectorized, and high-performance way.

Overview
--------

The function scans keyword arguments (``**kwargs``) for species-specific parameter groups
(e.g., ``lambda_1st``, ``f_2nd``, etc.), parses them by suffix, and constructs the full
parameter arrays for modeling. Only components with non-zero column densities are included.

Optimizations include:
- NumPy vectorization
- Reduced memory overhead
- Filtering out unused (zero-column-density) components

Accepted Parameter Groups
-------------------------

The following parameters must be provided for each species, with a suffix (e.g., `1st`, `2nd`, etc.):

- ``lambda_<suffix>``: Line center wavelengths
- ``f_<suffix>``: Oscillator strengths
- ``gamma_<suffix>``: Natural broadening (damping constants)
- ``b_<suffix>``: Doppler parameters (one per component)
- ``N_<suffix>``: Column densities (one per component)
- ``v_rad_<suffix>``: Radial velocities (km/s)

Each group must be self-consistent in length and structure.

Returns
-------

The output is passed directly from ``mother_function`` and typically represents an absorption
model evaluated over the provided wavelength grid.

Raises
------

- ``ValueError`` if any parameter group is missing or misaligned.
