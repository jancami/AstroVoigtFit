import numpy as np
from scipy.special import wofz
def voigt_profile(x, sigma, gamma):
    """
    Function to return the value of a (normalized) Voigt profile centered at x=0
    and with (Gaussian) width sigma and Lorentz damping (=HWHM) gamma.

    The Voigt profile is computed using the scipy.special.wofz, which returns
    the value of the Faddeeva function.


    WARNING
    scipy.special.wofz is not compaible with np.float128 type parameters.

    Args:
        x (float64): Scalar or array of x-values
        sigma (float64): Gaussian sigma component
        gamma (float64): Lorentzian gamma (=HWHM) component

    Returns:
        ndarray: Flux array for given input

    """

    z = (x + 1j * gamma) / sigma / np.sqrt(2)

    return np.real(wofz(z)) / sigma / np.sqrt(2 * np.pi)
