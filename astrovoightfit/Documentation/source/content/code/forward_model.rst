Forward Model
==============

Overview
--------

The `mother_function` models absorption lines using Voigt profiles and simulates a normalized absorption spectrum over a given wavelength grid. It takes into account Doppler and Lorentzian broadening, column density, and instrumental resolution, then returns a smoothed, interpolated model spectrum.



Function Definition
-------------------

.. code-block:: python

    def mother_function(wavegrid, lambda0=0.0, f=0.0, gamma=0.0, b=0.0, 
                        N=0.0, v_rad=0.0, v_resolution=0.0, map=map, n_step=25,
                        composition=None, debug=False):

        # Pre-compute all Voigt FWHMs at once
        Voigt_FWHM = VoigtFWHM(lambda0, gamma, b)
        FWHM2use = np.min(np.append(Voigt_FWHM, v_resolution))

        # Optimize grid creation
        xgrid_test = np.asarray(wavegrid)
        dv_xgrid = np.median(np.diff(xgrid_test)) / np.mean(xgrid_test) * cst.c.to("km/s").value
        n_step = max(7, np.ceil(FWHM2use / dv_xgrid))

        # Vectorized wavelength range calculation
        v_stepsize = FWHM2use / n_step
        bluewaves = lambda0 * (1.0 + (v_rad - 8.5 * Voigt_FWHM) / cst.c.to("km/s").value)
        redwaves = lambda0 * (1.0 + (v_rad + 8.5 * Voigt_FWHM) / cst.c.to("km/s").value)

        minwave = max(np.min(bluewaves), wavegrid.min())
        maxwave = min(np.max(redwaves), wavegrid.max())

        # Create reference grid
        n_v = int(np.ceil((maxwave - minwave) / minwave * cst.c.to("km/s").value / v_stepsize))
        refgrid = minwave * (1.0 + np.arange(n_v) * v_stepsize / cst.c.to("km/s").value)

        # Process all lines at once if possible
        allcomponents = np.zeros((n_v, len(lambda0)))

        for lineloop in range(len(lambda0)):
            dv = getVGrid(lambda0[lineloop], gamma[lineloop], b[lineloop], v_resolution, n_step)
            thiswavegrid = lambda0[lineloop] * (1.0 + dv / cst.c.to("km/s").value)
            tau = voigt_optical_depth(
                thiswavegrid,
                lambda0=lambda0[lineloop],
                b=b[lineloop],
                N=N[lineloop],
                f=f[lineloop],
                gamma=gamma[lineloop],
                v_rad=v_rad[lineloop],
            )
            vel = dv + v_rad[lineloop]
            thiswavegrid = lambda0[lineloop] * (1.0 + vel / cst.c.to("km/s").value)
            interpolationfunction = interp1d(thiswavegrid, tau, kind="cubic", fill_value="extrapolate")
            tau_grid = interpolationfunction(refgrid)
            tau_grid[np.where(refgrid > np.max(thiswavegrid))] = 0
            tau_grid[np.where(refgrid < np.min(thiswavegrid))] = 0
            allcomponents[:, lineloop] = tau_grid

        # Combine all lines
        tau = np.sum(allcomponents, axis=1)
        AbsorptionLine = np.exp(-tau)

        # Apply instrumental smoothing
        smooth_sigma = fwhm2sigma(v_resolution) / v_stepsize
        gauss_smooth = gaussian_filter(AbsorptionLine, sigma=smooth_sigma)

        # Interpolate final result
        interpolated_model = np.interp(wavegrid, refgrid, gauss_smooth, left=1, right=1)

        return interpolated_model

Explanation of Each Step
------------------------

**1. Voigt FWHM Calculation**

- Computes the full width at half maximum of the Voigt profile for each line:
  - Includes Doppler and Lorentzian broadening.
  - Uses the smallest FWHM to define step size.

**2. Velocity Step and Oversampling**

- `dv_xgrid`: Determines the average velocity spacing in `wavegrid`.
- `n_step`: Ensures there are enough samples per FWHM to avoid undersampling.

**3. Wavelength Range and Grid Setup**

- Calculates blue and red shifted limits around each line to find the valid wavelength range.
- Constructs a shared high-resolution reference grid `refgrid` for all lines.

**4. Line-by-Line Voigt Profile Calculation**

- For each line:
  - Computes velocity grid `dv`
  - Calculates optical depth `tau(λ)`
  - Converts velocity grid to wavelength
  - Interpolates optical depth to `refgrid`

**5. Optical Depth Summation**

- All interpolated τ arrays are summed along axis 1 to obtain the total optical depth.

**6. Radiative Transfer**

- Converts optical depth to transmission using:
  - `F(λ) = exp(-τ(λ))`

**7. Instrumental Smoothing**

- Uses a Gaussian filter to simulate instrumental broadening with a width defined by `v_resolution`.

**8. Final Interpolation**

- Uses fast 1D linear interpolation (`np.interp`) to match the original input `wavegrid`.

Why This Version is Faster
--------------------------

This optimized implementation is faster than previous versions due to the following changes:

**1. Vectorized Input Handling**

.. code-block:: python

    Voigt_FWHM = VoigtFWHM(lambda0, gamma, b)

Handles arrays from the start, avoiding repeated reshaping.

**2. Efficient Grid Construction**

.. code-block:: python

    minwave = max(np.min(bluewaves), wavegrid.min())
    maxwave = min(np.max(redwaves), wavegrid.max())

Reduces the wavelength grid size by focusing only on the relevant range.

**3. Fast Final Interpolation**

.. code-block:: python

    interpolated_model = np.interp(wavegrid, refgrid, gauss_smooth, left=1, right=1)

Uses `np.interp` instead of slower cubic interpolation.

**4. Preallocation of Arrays**

.. code-block:: python

    allcomponents = np.zeros((n_v, len(lambda0)))

Avoids repeated memory allocations inside the loop.

**5. Minimal Branching / Debug Logic**

All debug printing is commented out or optional, streamlining execution.

Conclusion
----------

This function is a performance-optimized implementation of Voigt-based absorption line modeling, using physics-based profiles and fast numerical schemes. It is especially useful in iterative fitting routines where model evaluation speed is critical.
