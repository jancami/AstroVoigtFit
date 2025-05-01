import numpy as np
import matplotlib.pyplot as plt
def errors_flux(waves, fluxes):
    """Takes the standard deviation of the first 50 normalized flux measurements 
        to estimate the error in the data set.

        Args:
            waves(ndarray): The wavelength data for the sightline
            fluxes(ndarray): The normalized flux for the sightline

        Returns:
            errors: ndarray
                Array of the estimated errors from the first 50 data points"""
    
    #standard deviation of the first 50 fluxes
    errors = np.std(fluxes[:50])
    #plot what the first fifty points look like, to see if its a somewhat random distribution
    plt.plot(waves[:50],fluxes[:50])
    return errors


