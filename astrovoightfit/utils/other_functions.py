import numpy as np
import astropy.constants as cst
from scipy.special import wofz
from scipy.interpolate import interp1d
from scipy.ndimage import gaussian_filter


def fwhm2sigma(fwhm):
    """
    Simple function to convert a Gaussian FWHM to Gaussian sigma.
    """
    sigma = fwhm / (2.0 * np.sqrt(2.0 * np.log(2.0)))
    return sigma


def VoigtFWHM(lambda0, gamma, b):
    """
    Calculate FWHM of Voigt using the formula by Olivero.
    See the bottom of wiki Voigt Profile page

    :param lambda0: center wavelength in AA
    :param gamma: Lorentzian gamma (=HWHM), in frequency
    :param b: Gaussian b in km/s.
    :return: Gau_FWHM, Lor_FWHM, V_FWHM, all in km/s
    """

    # gamma_AA = gamma * lambda0**2 / cst.c.to("angstrom/s").value
    gamma_kms = gamma * lambda0 * 1e-13
    Lor_FWHM = gamma_kms * 2

    # FWHM = 2*sigma*sqrt(2*ln2)
    Gau_FWHM = 2.35482 * b

    # fV = fL/2 + sqrt(fL^2/4 + fG^2)
    Voigt_FWHM = 0.5 * Lor_FWHM + np.sqrt(Lor_FWHM ** 2 / 4 + Gau_FWHM ** 2)

    return Voigt_FWHM


def getVGrid(lambda0, gamma, b, v_resolution, n_step):
    """
    Calculate v grid to be used for Voigt profile
    Size of grid is determined by the greater of FWHM_Voigt and v_resolution

    :return: dv, velocity grid
    """

    Voigt_FWHM = VoigtFWHM(lambda0, gamma, b)
    FWHM2use = np.max([Voigt_FWHM, v_resolution])
    v_stepsize = FWHM2use / n_step
    dv = np.arange(
        start=-8.5 * FWHM2use, stop=8.5 * FWHM2use, step=v_stepsize
    )
    return dv


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


def voigt_optical_depth(wave, lambda0=0.0, b=0.0, N=0.0, f=0.0, gamma=0.0, v_rad=0.0):
    """
    Function to return the value of a Voigt optical depth profile at a given wavelength, for a line
    centered at lambda0.

    Args:
        wave (float64): Wavelength at which optical depth is to be calculatd (in Angstrom)
        lambda0 (float64): Central (rest) wavelength for the absorption line, in Angstrom.
        b (float64): The b parameter (Gaussian width), in km/s.
        N (float64): The column density (in cm^{-2})
        f (float64): The oscillator strength (dimensionless)
        gamma (float64): Lorentzian gamma (=HWHM) component
        v_rad (float64): Radial velocity of absorption line (in km/s)

    Returns:
        float64: Optical Depth at wave.

    """

    # All we have to do is proper conversions so that we feed the right numbers into the call
    # to the VoigtProfile -- see documentation for details.
    nu = cst.c.to("angstrom/s").value / wave
    nu0 = cst.c.to("angstrom/s").value / lambda0
    sigma = (b * 1e13) / lambda0 / np.sqrt(2)
    gamma_voigt = gamma / 4 / np.pi
    tau_factor = (N * np.pi * cst.e.esu ** 2 / cst.m_e.cgs / cst.c.cgs * f).value

    # print("Nu0 is:        " + "{:e}".format(nu0))
    # print("Sigma is:      " + "{:e}".format(sigma))
    # print("Gamma is:      " + "{:e}".format(gamma_voigt))
    # print("Tau_factor is: " + "{:e}".format(tau_factor))

    # Transform this into a frequency grid centered around nu0

    ThisVoigtProfile = voigt_profile(nu - nu0, sigma, gamma_voigt)
    tau = tau_factor * ThisVoigtProfile

    return tau