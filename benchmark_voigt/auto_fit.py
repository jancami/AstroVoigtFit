import lmfit.model

import edibles.utils.benchmark_voigt.errors as e
import edibles.utils.benchmark_voigt.fits as fi
import edibles.utils.benchmark_voigt.initial_guesses as ig
import edibles.utils.benchmark_voigt.populate_params as pp
import edibles.utils.benchmark_voigt.read_data as rd
import matplotlib.pyplot as plt
from edibles.utils import voigt_profile as vp
import numpy as np

#get data for sightlines
resolution, star_name, central_wavelength = rd.read_benchmark_data(4226.728,4226)
f = 1.77
g = 2.200e+8
nstep = 25

for i in range(len(star_name)):
    #get welty and hobbs data
    b_WH,N_WH,v_rad_WH = rd.read_WH_data(star_name[i])
    #get star data
    wavelengths, normflux = rd.get_lamb_flux(star_name[i],4226)
    #initial parameter guesses
    b, N, v_rad = ig.initial_parm_guess(wavelengths,normflux, central_wavelength,resolution)
    #errors
    errors = e.errors_flux(wavelengths,normflux)
    #all of the components
    fitcomp,b,N,v_rad,delta_bic = fi.our_fit(wavelengths,normflux,central_wavelength,f,g,b,N, v_rad, resolution[i],nstep,errors, star_name[i],b_WH,N_WH)
    whfit = fi.welty_hobbs(wavelengths, normflux, central_wavelength, f, g, b_WH, N_WH, v_rad_WH, resolution[i], nstep,
                           errors)
    fit = pp.populate_params(fitcomp,b,delta_bic,b_WH,N_WH,v_rad_WH, wavelengths, normflux,central_wavelength,f,g,resolution[i],nstep,errors,star_name[i],whfit,errors)


