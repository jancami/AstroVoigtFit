from functions import PYTHONDIR
from pathlib import Path
import pandas as pd
import numpy as np
from edibles.utils.benchmark_voigt.parameter_modelling import file_reader


def read_benchmark_data(wavelength, sightline):
    """This function reads the csv files relating to each sightline, by element/molecule
        Args:
            wavelength(float64): the wavelengths in angstroms from the mail.txt file
            sightline(int): the rounded version of each wavelength ie. 5895
        Returns:
            ndarrays: containing the resolution,starname and central wavelength from each file"""
    
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
    """reads the welty and hobbs data
        Args:
            starname(string): the star name produced from read_benchmark_data
        Returns:
            b_WH: 1dnarray 
            an array containing all the broadening parameters for each component
            N_WH: 1dnarray 
            an array containing all the column densities for each component
            v_rad_WH: 1dnarray 
            an array containing all the radial velocities for each component"""
    
    # define path to file which contains the names of each star and the resolution of intraments used in survey
    folder = Path(PYTHONDIR + "/data")
    filename = folder / "voigt_benchmarkdata" / 'parameter_modelling_data' / (starname + '.txt')

    #define headers the data is under
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
    """ Will read in any file from the voigt_benchmark folder once given the name of the file.
        Args:
            starname(string): the star name produced from read_benchmark_data
            sightline(int): the wavelength used in read_benchmark_data

        Returns:
            wavelengths: 1dnarray 
            array containing all of the measured wavelengths for a given file
            normflux: 1dnarray 
            array containing all of the normalized flux values for a given file"""
    
    #casting the data from the file being read into arrays
    wavelengths, normflux = np.asarray(file_reader(starname + '.'+str(sightline) + '.txt'))
    
    return wavelengths, normflux





