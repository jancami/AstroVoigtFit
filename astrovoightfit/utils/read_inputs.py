from main_run import astrovoigtfit_run, continuum_fit
import matplotlib.pyplot as plt

#user inputs
star = "HD 183143"

#interested molecules
molecule = ['13CH+_4032','12CH+_4032']  # List of molecules to analyze

# file number of the observations
file_no = 0

# wavelength range for the fit  *** it's a highly sensitive parameter, so choose wisely
wave_range = [4231.5, 4233.5]


# absorption range is the region where continuum is not present i.e the region where the absorption lines are present
absorption_range = (4232.05, 4232.8)


# guessed species parameters for the molecules to fit the spectra 
species_params = {
    0: {  # 13CH+ 
        'b': [2, 2],
        'N': [1e13/70, 1e13/70],
        'v_rad': [-11, 4]
    },
    1: {  # 12CH+ 
        'b': [2, 2],
        'N': [1e13, 1e13],
        'v_rad': [-11, 4]
    }
}


# before running lets check if the continuum fitted properly or not 
wave,continuum_normalized_flux ,flux, continuum = continuum_fit(star, molecule,file_no, wave_range, absorption_range)

plt.figure(figsize=(8, 8))
plt.subplot(2, 1, 1)
plt.plot(wave, flux, 'b-', label='Original Flux')
plt.plot(wave, continuum, 'r-', label='Continuum Fit')
plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3, label='Absorption Region')
plt.xlabel('Wave')
plt.ylabel('Flux')
plt.legend()
plt.title('Continuum Fitting')

# plotting continuu, Normalized flux
plt.subplot(2, 1, 2)
plt.plot(wave, continuum_normalized_flux, 'g-', label='Normalized Flux')
plt.axhline(1.0, color='k', linestyle='--', label='Continuum = 1.0')
plt.axvspan(absorption_range[0], absorption_range[1], color='gray', alpha=0.3)
plt.xlabel('Wavelength')
plt.ylabel('Normalized Flux')
plt.legend()
plt.title('Normalized Spectrum')
plt.show()


user_analysis = input("Do you want to run the fitting model? (yes/no): ").strip().lower()
    
if user_analysis == 'yes' or user_analysis == 'y':
    print("running the fitting model")
    astrovoigtfit_run(star,molecule, wave_range,species_params,absorption_range,file_no)
else: 
    print("Please rerun using your desired parameters.")