from setuptools import setup, find_packages

setup(
    name="AstroVoigtFit",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy",
        "matplotlib",
        "specutils"
        "astropy>=4.0",
        "pandas",
        "lmfit",
    ],
)