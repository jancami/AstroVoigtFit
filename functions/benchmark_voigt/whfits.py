import numpy as np
import matplotlib.pyplot as plt
from functions.benchmark_voigt.voigt_profile import *
from functions.benchmark_voigt import parameter_modelling as pm


import functions.benchmark_voigt.errors as e
import functions.benchmark_voigt.read_data as rd
from functions.benchmark_voigt.benchmark_fit import *


#this file produces just the welty and hobbs data, similar to the way auto_fit.py does it
#just need to switch the variables in the calls of read_benchmark_data and get_lamb_flux, as well as f and g according to the wavelengths
#in the mail.txt file in the parameter_modelling_data folder
#could automate this variable changing

#get data for sightlines
resolution, star_name, central_wavelength = rd.read_benchmark_data(5895.9242,5895)
f = 3.270e-1
g =6.280e+7
nstep = 25

for i in range(len(star_name)):
    #get welty and hobbs data
    b_WH,N_WH,v_rad_WH = rd.read_WH_data(star_name[i])
    #get star data
    wavelengths, normflux = rd.get_lamb_flux(star_name[i],5895)

    #errors
    errors = e.errors_flux(wavelengths,normflux)



    #calling fit_multi....
    fittedmodelman = fit_multi_voigt_absorptionlines(wavegrid = wavelengths, ydata = normflux, restwave = central_wavelength, f = f,
                                     gamma = g, b = b_WH, N = N_WH, v_rad = v_rad_WH, v_resolution = 0.50, n_step=nstep, std_dev =errors)
    
    #this step prints out each individual component used in the fit
    for comp in range(len(b_WH)):
        # print(comp)
        indivfit = fit_multi_voigt_absorptionlines(wavegrid=wavelengths,
                                                      ydata=normflux,
                                                      restwave=central_wavelength,
                                                      f=f,
                                                      gamma=g,
                                                      b=b_WH[comp],
                                                      N=N_WH[comp],
                                                      v_rad=v_rad_WH[comp],
                                                      v_resolution=0.50,
                                                      n_step=nstep,
                                                      std_dev=errors)
        fig = plt.subplot(111)
        plt.plot(wavelengths, indivfit.best_fit, label=comp)
        plt.legend()
        plt.show()

    #plot statement
    plt.plot(wavelengths, normflux, label='data')
    plt.plot(wavelengths,fittedmodelman.best_fit, label = 'fitted model')
    #the statistical values
    print(fittedmodelman.redchi)
    print(fittedmodelman.bic)
    print(fittedmodelman.aic)
    plt.legend()
    plt.show()
