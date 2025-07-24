   ___         _              __      ___     _ _ _ _ _       
  / _ \___ ___| |_ ___  _ __  \ \    / (_) __| (_) | |_ _  _ 
 | | | / _ (_-<  _/ _ \| '_ \  \ \/\/ /| |/ _` | | | |  _| || |
 | |_| \_,__/__/\__\___/| .__/   \_/\_/ |_|\__,_|_|_|_|\__|\_, |
  \___/                  |_|                              |__/  

AstroVoigtFit
=============

A Python package for fitting astrophysical absorption spectra using Voigt profiles.

Key Features:
-------------

- `read_inputs.py`: The **master function** to run Voigt profile fits from input files.
- `continuum.py`: Fits the **continuum** of the spectrum before line fitting.
- `astrovoigtfit.py`: Core module to **fit individual species** using Voigt profiles.

Usage:
------

```bash
python astrovoightfit/utils/read_inputs.py

