import numpy as np
import matplotlib.pyplot as plt
from functions.benchmark_voigt import voigt_profile as vp
from functions.benchmark_voigt.parameter_modelling import file_reader
import astropy.constants as cst
import scipy.signal as ss
from edibles import PYTHONDIR
from pathlib import Path
import pandas as pd

def find_nearest(array, value):
    """
    finds the value in an array which is closest to the desired value
    paramaters :
        array : list/array/numpy array
            an array of vaules which we plan to extrcat a single values from
        value : int
            the value which we wish to be close to
    return:
        idx : int
            the position in the array of the nearest value
    """
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx

def find_FWHM (x_data, y_data):
    """
    Will calculate the full width half maximum from data entered
    parametres:
        x_data: numpy array
        y_data: numpy array
    return:
        fwhm: float
            value of the fwhm in x_data units
        peak: interger
            position of peak value within y_data
    """
    # determine the index of the peak in y_data
    peak = np.argwhere(y_data == np.min(y_data))[0,0]
    # Find half of the peak value
    half_peak = y_data[peak] + (1 - y_data[peak])*0.5
    #find nearst data points to half_peak on both sides of the peak
    left_side = find_nearest(y_data[:peak], half_peak)
    right_side = find_nearest(y_data[peak:], half_peak) + len(y_data[:peak])
    # calculae the full width half maximum (fwhm)
    fwhm = x_data[right_side] - x_data[left_side]
    return fwhm, peak

