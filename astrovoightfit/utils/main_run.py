from edibles.utils.edibles_oracle import EdiblesOracle
from edibles.utils.edibles_spectrum import EdiblesSpectrum
from astrovoigtfit import *

# calling standard python libraries
import numpy as np
import matplotlib.pyplot as plt
import re


def observations(star,molecule,file_no, wave_range,species_file='species.txt'):
    
    with open(species_file) as f:
        lines = f.readlines()

    headers = lines[0].split()
    line_index = headers.index('line')

    for line in lines[1:]:
        parts = line.split()
        if parts[0] == molecule:
            line_val = parts[line_index].strip('[]')
            line_val = float(line_val)
            # print(line_val)
            # print(wave_range)
            
    pythia = EdiblesOracle()
    List = pythia.getFilteredObsList(object=[star], MergedOnly=False, Wave=line_val)

    # print(List)
    test = List.values.tolist()
    # print(test)
    
    filename = test[file_no]
    print(filename)
    wrange = wave_range
    sp = EdiblesSpectrum(filename)
    sp.getSpectrum(wrange[0],wrange[1])
    wave = sp.bary_wave
    flux = sp.bary_flux
    idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
    wave = wave[idx]
    flux = flux[idx]
    flux = flux / np.median(flux)

    # Filter the spectrum within the specified wavelength range
    idx = np.where((wave > wrange[0]) & (wave < wrange[1]))
    wave = wave[idx]
    flux = flux[idx]

    # Normalize the flux
    flux = flux / np.median(flux)
    
    return wave, flux
    
    
def continuum_fit(star, molecule,file_no, wave_range, absorption_range):
    
    # Get the observed spectrum
    wave, flux = observations(star, molecule[0], file_no, wave_range)

    absorption_range = absorption_range  # Wavelength range of the absorption feature
    degree = 3  # Degree of Chebyshev polynomial

    #  Fit continuum and normalize
    continuum_normalized_flux, continuum, poly, std_dev = fit_continuum(
        wave, flux, absorption_range, degree, return_std=True
    )  
    
    return wave,continuum_normalized_flux ,flux, continuum  

    # # plotting the normalized spectrum and the continuum curve
    # plt.figure(figsize=(8, 8))
    # plt.subplot(2, 1, 1)
    # plt.plot(wave, flux, 'b-', label='Original Flux')
    # plt.plot(wave, continuum, 'r-', label='Continuum Fit')
    # plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3, label='Absorption Region')
    # plt.xlabel('Wave')
    # plt.ylabel('Flux')
    # plt.legend()
    # plt.title('Continuum Fitting')

    # # fitting continuum Normalized flux
    # plt.subplot(2, 1, 2)
    # plt.plot(wave, continuum_normalized_flux, 'g-', label='Normalized Flux')
    # plt.axhline(1.0, color='k', linestyle='--', label='Continuum = 1.0')
    # plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3)
    # plt.xlabel('Wavelength')
    # plt.ylabel('Normalized Flux')
    # plt.legend()
    # plt.title('Normalized Spectrum')


def _parse_bracketed_list(s):
    inner = s.strip().strip('[]')
    if not inner:
        return []
    return [float(x.strip()) for x in inner.split(',') if x.strip()]

def get_species_params(species_file, species_params, molecule):
    with open(species_file) as f:
        lines = f.readlines()

    for idx, species in enumerate(molecule):
        for line in lines[1:]:
            if not line.strip():
                continue
            # species name is first whitespace-separated token
            parts = line.split()
            if parts[0] != species:
                continue

            # extract the three bracketed columns (lambda, f_value, gamma)
            bracketed = re.findall(r'\[.*?\]', line)
            if len(bracketed) < 3:
                raise ValueError(f"Line for species {species} does not have expected bracketed fields: {line!r}")

            lambda_vals = _parse_bracketed_list(bracketed[0])
            f_vals = _parse_bracketed_list(bracketed[1])
            gamma_vals = _parse_bracketed_list(bracketed[2])

            print(f"Species: {species}, Lambda: {lambda_vals}, f: {f_vals}, Gamma: {gamma_vals}")

            reordered = {
                'lambda': lambda_vals,
                'f': f_vals,
                'gamma': gamma_vals,
            }

            for key, value in species_params[idx].items():
                reordered[key] = value

            species_params[idx] = reordered
            break  # found species line, move to next

    return species_params




def astrovoigtfit_run(star,molecule, wave_range,species_params,absorption_range,file_no,species_file='species.txt'):
    
    # Get the observed spectrum and fit the continuum
    
    wave, flux = observations(star,molecule[0],file_no, wave_range,species_file)
    species_param_updated = get_species_params(species_file, species_params, molecule)
    # print("Species parameters:", species_param_updated)
    
    absorption_range = absorption_range  # Wavelength range of the absorption feature
    degree = 3  # Degree of Chebyshev polynomial

    #  Fit continuum and normalize
    continuum_normalized_flux, continuum, poly, std_dev = fit_continuum(
        wave, flux, absorption_range, degree, return_std=True
    )     
    

    # fitting the data using astro_voigt_fit function
    fitresult= astro_voigt_fit(
        wavegrid=wave, 
        ydata=continuum_normalized_flux,
        species_params=species_param_updated,
        v_resolution=3, 
        n_step=25, 
        std_dev=0.0014
    )
    fitresult.params.pretty_print() #printing the fitting parameters



    print("chi-square value ",fitresult.chisqr)
    print("reduced chi-square value ",fitresult.redchi)
    print("FITTING RESULT :", fitresult.success)


    plt.plot(wave,fitresult.best_fit,color ='purple',label ="fit")
    plt.plot(wave, continuum_normalized_flux,color ='gray',label ='data',alpha = 0.7)
    plt.xlabel("Wavelength ($\AA$)")
    plt.ylabel("Normalised flux")
    plt.title("Multi cloud single line model fit for CH+",color = 'darkgreen')
    plt.grid()
    plt.legend()
    plt.show()