Other Functions
===============

This module provides utility functions for working with Voigt profiles and Gaussian–Lorentzian line shapes. These are used internally by the `mother_function` to compute velocity grids, FWHMs, and optical depth profiles.

Function List
-------------

.. autofunction:: other_functions.fwhm2sigma
.. autofunction:: other_functions.VoigtFWHM
.. autofunction:: other_functions.getVGrid
.. autofunction:: other_functions.voigt_profile
.. autofunction:: other_functions.voigt_optical_depth


Function Details
----------------

fwhm2sigma
~~~~~~~~~~

.. code-block:: python

    fwhm2sigma(fwhm)

Convert a Gaussian full width at half maximum (FWHM) into the Gaussian standard deviation (sigma).

**Parameters:**

- `fwhm` (float): Full width at half maximum.

**Returns:**

- `sigma` (float): Standard deviation of the Gaussian profile.


VoigtFWHM
~~~~~~~~~

.. code-block:: python

    VoigtFWHM(lambda0, gamma, b)

Compute the Voigt profile's FWHM using Olivero's approximation.

**Parameters:**

- `lambda0` (float or array): Central wavelength(s) in Angstrom.
- `gamma` (float): Lorentzian width (HWHM).
- `b` (float): Gaussian width in km/s.

**Returns:**

- `Voigt_FWHM` (float): Voigt profile FWHM in km/s.


getVGrid
~~~~~~~~

.. code-block:: python

    getVGrid(lambda0, gamma, b, v_resolution, n_step)

Compute the velocity grid for evaluating the Voigt profile, based on either instrumental resolution or the natural FWHM of the line.

**Returns:**

- `dv` (array): Velocity grid in km/s.


voigt_profile
~~~~~~~~~~~~~

.. code-block:: python

    voigt_profile(x, sigma, gamma)

Compute the normalized Voigt profile using the Faddeeva function (via `scipy.special.wofz`).

**Returns:**

- `flux` (array): Normalized profile values.


voigt_optical_depth
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

    voigt_optical_depth(wave, lambda0, b, N, f, gamma, v_rad)

Compute the Voigt-profile-based optical depth for a given wavelength array and line parameters.

**Parameters:**

- `wave` (array): Wavelength array (in Angstrom).
- `lambda0` (float): Central wavelength.
- `b` (float): Doppler broadening (Gaussian width) in km/s.
- `N` (float): Column density (in cm⁻²).
- `f` (float): Oscillator strength.
- `gamma` (float): Lorentzian width (HWHM).
- `v_rad` (float): Radial velocity shift (in km/s).

**Returns:**

- `tau` (array): Optical depth at each input wavelength.


Notes
-----

- `voigt_profile` internally uses `scipy.special.wofz`, which may not support `np.float128`. Stick to `np.float64` inputs.
- The `VoigtFWHM` function uses an empirical approximation valid for most astrophysical conditions.
