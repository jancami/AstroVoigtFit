from functions.benchmark_voigt import voigt_profile as vp
import scipy.signal as ss
import numpy as np
import astropy.constants as cst
import matplotlib.pyplot as plt
import functions.benchmark_voigt.initial_guesses as ig
from functions.benchmark_voigt import auto_voigt_fit as av




def welty_hobbs(wavelen,nflux,cenw,f,g,b,n,v,vres,nstep,stddev):
    """ calls the fit_multi_voigt_absorptionlines function with the welty and hobbs parameters.
        Args:
            wavelen(ndarray): wavelength array (angstroms)
            nflux(ndarray): normalized flux array
            cenw(float64): central wavelength (angstroms)
            f(float64): oscilator strength
            g(float64): Lorentzian gamma
            b(float64): the b parameter (gaussian width) (km/s)
            n(float64):  The column density (in cm^{-2})
            v (float64): Radial velocity of absorption line (in km/s)
            vres(float64): Instrument resolution in velocity space (in km/s)
            nstep (int): no. of point per FWHM length, governing sampling rate and efficiency
            stddev: The errors used for the fit_multi_voigt_absorptionlines call

        Returns:
            whfit: an lmfit model instance"""
    
    whfit = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelen,
                                               ydata=nflux,
                                               restwave=cenw,
                                               f=f,
                                               gamma=g,
                                               b=b,
                                               N=n,
                                               v_rad=v,
                                               v_resolution=vres,
                                               n_step=nstep,
                                               std_dev=stddev)
    return whfit


def our_fit(wavelen,nflux,cenw,f,g,bour,nour,vour,vres,nstep,stddev,starname,bwh,nwh):
    """ calls the fit_multi_voigt_absorptionlines function with our inital parameter guesses from inital_guesses.py. 
        Also determines whether to keep the fit going based off of the reduced chi squared, AIC and BIC.
        
        Args:
            wavelen(ndarray): wavelength array (angstroms)
            nflux(ndarray): normalized flux array
            cenw(float64): central wavelength (angstroms)
            f(float64): oscilator strength
            g(float64): Lorentzian gamma
            bour(float64): the b parameter (gaussian width) (km/s)
            nour(float64):  The column density (in cm^{-2})
            vour (float64): Radial velocity of absorption line (in km/s)
            vres(float64): Instrument resolution in velocity space (in km/s)
            nstep (int): no. of point per FWHM length, governing sampling rate and efficiency
            stddev: The errors used for the fit_multi_voigt_absorptionlines call
            starname(string): Name of the star to label the plots
            bwh(ndarry): the b parameter from welty and hobbs
            nwh(ndarray): the N parameter fromm welty and hobbs

        Returns:
            fit: an lmfit model instance
            b: ndarray
            the b parameter array produced by the fitting algorthm
            n: ndarray
            the N parameter array produced by the fitting algorthm
            v_rad: ndarray
            the v_rad parameter array produced by the fitting algorthm
            delta_bic: the change in BIC between the final fit, and its previous fit
            """
    
    fit1 = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelen,
                                             ydata=nflux,
                                             restwave=cenw,
                                             f=f,
                                             gamma=g,
                                             b=bour,
                                             N=nour,
                                             v_rad=vour,
                                             v_resolution=vres,
                                             n_step=nstep,
                                             std_dev=stddev)


    # min_continuum, max_contiumm = np.argwhere()
    residual_array = nflux - fit1.best_fit

    # find the peak position of the residuals plot
    peak_position = wavelen[np.argwhere(np.abs(residual_array) == np.max(np.abs(residual_array)))[0, 0]]

    # create an array of all minima in the residuals plot
    mins = wavelen[ss.argrelextrema(residual_array, np.less)]

    print('peak ', peak_position)
    # print('mins', find_nearest(mins[np.argwhere(mins > peak_position)], peak_position))

    # find where the 2 nearest local minimum to the peak is
    lower_near_lambda = mins[av.find_nearest(mins[np.argwhere(mins < peak_position)], peak_position)]
    upper_near_lambda = mins[
        av.find_nearest(mins[np.argwhere(mins > peak_position)], peak_position) + np.argwhere(mins > peak_position)[0, 0]]

    # find the fwhm and peak positon of the 1 component fit
    fwhm, peak_position = av.find_FWHM(wavelen, fit1.best_fit)
    # limit the data being analysed to within 3 fwhm of the central peak
    #ight_continuum = av.find_nearest(wavelen, wavelen[peak_position] + (5 * fwhm))
    #left_continuum = av.find_nearest(wavelen, wavelen[peak_position] - (5 * fwhm))
    #wavelen = wavelen[left_continuum:right_continuum]
    #nflux = nflux[left_continuum:right_continuum]
    # fit.best_fit = fit.best_fit[left_continuum:right_continuum]
    print('reduced chi-square: ', fit1.summary()['redchi'])
    # print('delta bic', delta_bic)
    # print('delta aic', delta_aic)
    # print('delta reduced chi squared:', delta_red_chi)

    # plotting results

    fig, (ax1, ax2, ax3,ax4) = plt.subplots(4, 1)
    # plt.xlim(7699,7699.7)
    ax1.set_xlim(4226.8, 4227)
    ax2.set_xlim(4226.8, 4227)
    ax3.set_xlim(4226.8, 4227)
    ax4.set_xlim(4226.8, 4227)
    ax4.plot(wavelen,fit1.init_fit)
    ax2.plot(wavelen, nflux, c='k', label='data')
    ax2.plot(wavelen, fit1.best_fit, 'b', label='1 component')
    ax1.plot(wavelen, nflux / fit1.best_fit, c='r', label='residuals')
    ax2.plot([lower_near_lambda] * 50, np.linspace(0, 2.2, 50), c='orange', linestyle='dashed', label='local minimas')
    ax2.plot([upper_near_lambda] * 50, np.linspace(0, 2.2, 50), c='orange', linestyle='dashed')
    #ax2.plot([peak_position] * 50, np.linspace(0, 2.2, 50), c='maroon', linestyle='dashed', label='peak residual')
    ax3.plot(wavelen, nflux - fit1.best_fit, label='actual residuals')
    ax1.legend(loc="upper right")
    ax2.legend(loc="upper right")
    ax3.legend(loc="upper right")
    # ax.plot(wavelengths[left_continuum:right_continuum], normflux[left_continuum:right_continuum], c = 'orange')
    fig.suptitle(starname)

    # saveing for use in project report
    # if star_name[i] == 'omiper.m95':
    # plt.savefig(
    # 'c:/Users/user/edibles/edibles/data/voigt_benchmarkdata/parameter_modelling_data/one_componet_fitting_for_document.png')
    plt.show()

    # assigning 'previous_bic' (bic = bayesian info criterium) this may be used in
    # component fitting in order to tell our program when we have a good fit
    previous_bic = fit1.summary()['bic']
    previous_aic = fit1.summary()['aic']

    # define reduced chi and delta_bic to allow us to enter the loop
    previous_red_chi = fit1.summary()['redchi']
    delta_bic = previous_bic
    delta_aic = previous_aic
    delta_red_chi = -1

    print("components chisqr BIC AIC")
    # minimum accepted previous_red_chi for working data
    # omiper.m95 previous_red_chi > 1
    # zetoph.lf previous_red_chi > 1.18
    # zetoph.k94 previous_red_chi > 1.47
    # need to write a constraint to stop the loop once it reaches the maximium number of parameters
    # count = len(b_WH) #1.47
    # (previous_red_chi > 1.47)
    fit = fit1

    v_rad = vour

    while (delta_bic < 0) | (delta_aic < 0) | (delta_red_chi < 0):#|(previous_red_chi > 6.5):
        # count -= 1
        # print(count)
        # retrieve the previous reduced chi squared
        previous_red_chi = fit.summary()['redchi']

        # redifine previous fit before adding component to fit, previous fit will then be used outside of the while loop
        # and will be redifined as the best fit
        previous_fit = fit

        # array of residuals obtained by dividing the normflux by the flux from the last fit
        residual_array = nflux-fit.best_fit


        b,N,v_rad = ig.multicomp(wavelen,v_rad,residual_array,cenw,vres)

        # fit all of the components
        fit = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelen,
                                                 ydata=nflux,
                                                 restwave=cenw,
                                                 f=f,
                                                 gamma=g,
                                                 b=b,
                                                 N=N,
                                                 v_rad=v_rad,
                                                 v_resolution=vres,
                                                 n_step=nstep,
                                                 std_dev=stddev)
        # fit.params.pretty_print()
        indivfits = []
        for comp in range(len(b)):
        # print(comp)
            indivfit = vp.fit_multi_voigt_absorptionlines(wavegrid=wavelen,
                                                          ydata=nflux,
                                                          restwave=cenw,
                                                          f=f,
                                                          gamma=g,
                                                          b=b[comp],
                                                          N=N[comp],
                                                          v_rad=v_rad[comp],
                                                          v_resolution=vres,
                                                          n_step=nstep,
                                                          std_dev=stddev)
            indivfits.append(indivfit.best_fit)
            fig = plt.subplot(111)
            plt.plot(wavelen, indivfit.best_fit, label=comp)
            plt.plot(wavelen,(residual_array+1))
            plt.legend()
            plt.show()
        #print(indivfits)
        # calculate the difference in the bic from this fit and the previous fit, when this value is positive it means
        # that we have added too many components and should use the fit before this as th ideal fit
        delta_bic = fit.summary()['bic'] - previous_bic
        previous_bic = fit.summary()['bic']
        delta_red_chi = fit.summary()['redchi'] - previous_red_chi
        delta_aic = fit.summary()['aic'] - previous_aic
        previous_aic = fit.summary()['aic']

        # print('reduced chi-square: ',fit.summary()['redchi'])
        # print('bic', fit.summary()['bic'])
        # print('delta bic', fit.summary()['bic'] - previous_bic)

        # This is used for debugging just to check what is being fitted and when
        residual_array = nflux - fit.best_fit

        # find the peak position of the residuals plot
        # print('peak',wavelengths[np.argwhere(residual_array == np.max(residual_array))[0, 0]])
        # print('lambdas',(v_rad * central_wavelength / cst.c.to("km/s").value) + central_wavelength)
        print(len(v_rad), fit.summary()['redchi'], fit.summary()['bic'], fit.summary()['aic'])
        print(fit.summary()['chisqr'])
        print(np.sum((fit.best_fit - nflux) ** 2 / stddev ** 2))
        # print('delta bic', delta_bic)
        # print('delta aic', delta_aic)
        # print('delta reduced chi squared:', delta_red_chi)
        # print('previous reduced chi squared:', previous_red_chi)

        # print(fit.fit_report())

        # plotting data for each component, makes it easier to check the fitting process

        fig, (ax1, ax2, ax3,ax4) = plt.subplots(4, 1)
        # ax = fig.add_subplot(111)
        # ax1.set_ylim(-0.1,1.5)
        # ax2.set_ylim(-0.1, 1.5)
        # ax3.set_ylim(-0.1, 1.5)
        ax1.set_xlim(4226.8, 4227)
        ax2.set_xlim(4226.8, 4227)
        ax3.set_xlim(4226.8, 4227)
        ax4.plot(wavelen, fit.init_fit)
        ax2.errorbar(wavelen, nflux, yerr=np.array([stddev] * len(wavelen)), mfc='none', c='k',
                     label='Data', alpha=0.4)
        ax2.plot(wavelen, fit.best_fit, c='b', label=len(v_rad))
        ax1.plot(wavelen, nflux / fit.best_fit, label='residuals')
        ax2.plot(
            np.array([wavelen[np.argwhere(np.abs(residual_array) == np.max(np.abs(residual_array)))[0, 0]]] * 50),
            np.linspace(0, 2, 50), linestyle='dashed')
        ax3.plot(wavelen, nflux - fit.best_fit, label='actual residuals')
        for i in range(len(b)):
            ax2.plot(wavelen,indivfits[i])
        #ax2.plot(wavelen, indivfits[0], c = 'green')#test and fix/remove on wednesday
        #ax2.plot(wavelen, indivfits[1], c = 'red', alpha = 0)
        fig.suptitle(starname)
        ax1.legend(loc="upper right")
        ax2.legend(loc="upper right")
        ax3.legend(loc="upper right")
        plt.show()


    fit = previous_fit

    #fig = plt.figure(constrained_layout=True)
    #axs = fig.subplot_mosaic([['TopLeft', 'Right'], ['BottomLeft', 'Right']],
                             #gridspec_kw={'width_ratios': [2, 1]})
    #axs['TopLeft'].set_title('Component One')
    #axs['BottomLeft'].set_title('Component Two')
    #axs['Right'].set_title('Final Fit')
    #axs['TopLeft'].plot(wavelen, fit1.best_fit)
    #axs['TopLeft'].plot(wavelen, ((nflux - fit1.best_fit) + 1))
    #axs['TopLeft'].set_xlim(7699.0, 7699.6)
    #axs['BottomLeft'].plot(wavelen, indivfit.best_fit)
    #axs['BottomLeft'].plot(wavelen, ((nflux - fit.best_fit) + 1))
    #axs['BottomLeft'].set_xlim(7699.0, 7699.6)
    #axs['Right'].plot(wavelen, fit.best_fit)
    #axs['Right'].plot(wavelen, nflux, alpha=0.4, c='black')
    #axs['Right'].set_xlim(7699.0, 7699.6)
    #axs['TopLeft'].legend(loc="upper right")
    #plt.legend()
    #plt.show()

    # fig, (ax1, ax2) = plt.subplots(2, 1)
    # ax1.plot(wavelen,nflux, label = 'data', c= 'black', alpha = 0.4)
    # ax1.plot(wavelen,fit.best_fit, c = 'blue', label = 'Final Fit')
    # ax1.plot(wavelen, fit1.best_fit, c = 'green', label = 'First Component')
    # ax1.plot(wavelen, indivfit.best_fit, c = 'red', label = 'Second Component', alpha = 0.4)
    # ax2.plot(wavelen, ((nflux-fit1.best_fit)+1), c = 'green', label = 'First Component Residuals')
    # ax2.plot(wavelen, ((nflux - fit.best_fit) + 1), c = 'red', label = 'Second Component Residuals')
    # ax1.legend(loc = 'upper right')
    # ax2.legend(loc = 'upper right')
    # plt.legend()
    #plt.savefig("omiper.k94.png")
    # plt.show()

    # print(len(b))
    return fit,b,N,v_rad,delta_bic






