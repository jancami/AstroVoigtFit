import numpy as np
import matplotlib.pyplot as plt
def errors_flux(waves, fluxes):
    errors = np.std(fluxes[:50])
    plt.plot(waves[:50],fluxes[:50])
    return errors


