from functions.benchmark_voigt import voigt_profile as vp
import scipy.signal as ss
import numpy as np
import astropy.constants as cst
import matplotlib.pyplot as plt

from functions.benchmark_voigt import auto_voigt_fit as av

def initial_parm_guess(wave,flux,wavel):
    """this produces the inital guesses for b, N and v_rad for a one component fit.
        Args:
            wave: the wavelengths (angstroms)
            flux: the normalized flux
            wavel: the central wavelength (angstroms)
        Returns:
            b: ndarray
            array containing the initial b value
            N: ndarray
            array containing the initial N value
            v_rad: ndarray
            array containing the initial v_rad value"""
            
    b = np.array([0.25])
    N = np.array([1e12])

    # need good initial guess for v_rad
    # calc flux weighted wavelength and convert to v_rad
    wave_weight = sum((1 - flux) * wave / sum(1 - flux))
    print('wave weight', wave_weight)

    # convert the guess wave position into velocity using eq below
    v_rad_guess = (wave_weight - wavel) / wavel * cst.c.to("km/s").value
    print('v_rad', v_rad_guess)
    v_rad = np.array([v_rad_guess])

    return b,N,v_rad

def multicomp(wl,vr,res,cw,vres):
    """produces the initial guesses for b,N and v_rad for multiple component fits, 
    uses the absolute maximum of the residuals to determine v_rad.
        Args:
            wl: the wavelengths (angstroms)
            vr: radial velocity in km/s
            res: residuals
            cw: the central wavelength (angstroms)
            vres: instrumental resolution (km/s)
        Returns:
            b: ndarray
            array containing the initial b values
            N: ndarray
            array containing the initial N values
            v_rad: ndarray
            array containing the initial v_rad values"""
            
               
    # find the peak position of the residuals plot
    peak_position = wl[np.argwhere(np.abs(res) == np.max(np.abs(res)))[0, 0]]
    # peak_position = 1
    # create an array of all minima in the residuals plot
    mins = wl[ss.argrelextrema(res, np.less)]

    print('peak ', peak_position)
    # print('mins', find_nearest(mins[np.argwhere(mins > peak_position)], peak_position))

    # find where the 2 nearest local minimum to the peak is
    lower_near_lambda = mins[av.find_nearest(mins[np.argwhere(mins < peak_position)], peak_position)] - 0.05
    upper_near_lambda = mins[
                            av.find_nearest(mins[np.argwhere(mins > peak_position)], peak_position) +
                            np.argwhere(mins > peak_position)[
                                0, 0]] + 0.05

    # print('lower lambda', lower_near_lambda)
    # print('upper lambda', upper_near_lambda)

    # convert the previously fitted v_rads into wavelength for easy comparison with the upper and lower minimas
    min_lambdas = (vr * cw / cst.c.to("km/s").value) + cw

    # replace the v_rad nearest to the peak of the residuals data with the lower minima from above and then add
    # the upper minima to the  end of the array and sort the array into size order using .sort
    min_lambdas[av.find_nearest(min_lambdas, peak_position)] = lower_near_lambda
    min_lambdas = np.append(min_lambdas, upper_near_lambda)
    min_lambdas.sort()

    # convert all lambdas back to velocities for fitting
    v_rad = (min_lambdas - cw) / cw * cst.c.to("km/s").value + vres * 0.5

    # vdiff = [(v_rad[i]-v_rad[i-1]) for i in range(len(v_rad))]
    # v_rad += vdiff
    # make sure the b and N arrays are the right length
    b = np.array([0.85] * len(v_rad))
    N = np.array([7.5e10] * len(v_rad))

    # check that we are adding the right amount of components
    # print('components to be fitted', len(v_rad))

    #additional change to intial guesses after reaching 5 components as the components get much smaller
    #the fitting routine targets larger components first
    if len(v_rad) > 5:
        # v_rad[-1] = -5
        b[-1] = 0.45
        N[-1] = 1e10

    return b,N,v_rad
