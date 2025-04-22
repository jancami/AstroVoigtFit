from edibles import PYTHONDIR
from pathlib import Path
import pandas as pd
import numpy as np
from edibles.utils.benchmark_voigt.parameter_modelling import file_reader


def read_benchmark_data(wavelength, sightline):
    # define path to file which contains the names of each star and the resolution of intraments used in survey
    folder = Path(PYTHONDIR + "/data")
    filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / (str(sightline) +"files_for_auto_voigt.txt")

    # state what headers the desired data is under
    Headers = ["star_name", "resolution"]

    # read in the data
    file = pd.read_csv(
        filename,
        delim_whitespace=True,
        header=None,
        names=Headers,
        engine="python",
    )

    # collect the data from each header in easy to iterate form
    resolution = file['resolution']
    star_name = file["star_name"]
    central_wavelength = wavelength
    return resolution,star_name,central_wavelength

def read_WH_data(starname):
    # get parameters from Welty and Hobbs paper for comparison
    folder = Path(PYTHONDIR + "/data")
    filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / (starname + '.txt')

    Headers = ["b", "N", "v_rad"]
    star_parameters = pd.read_csv(filename,
                                  delim_whitespace=True,
                                  header=None,
                                  names=Headers,
                                  engine="python",
                                  )

    # collect the data from each header in easy to iterate form
    b_WH = np.asarray(star_parameters["b"])
    N_WH = np.asarray(star_parameters["N"])
    v_rad_WH = np.asarray(star_parameters["v_rad"])

    return b_WH,N_WH,v_rad_WH

def get_lamb_flux(starname, sightline):
    wavelengths, normflux = np.asarray(file_reader(starname + '.'+str(sightline) + '.txt'))
    #print(star_name)
    return wavelengths, normflux





