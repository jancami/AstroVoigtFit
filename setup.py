from setuptools import setup, find_packages

setup(
    name="AstroVoigtFit",
    version="0.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "scipy>=0.19",
        "matplotlib",
        "specutils"
        "astropy>=4.0.0",
        "pandas",
        "lmfit",
        "cycler",
        "kiwisolver",
        "pyparsing",
        "python-dateutil",
        "pytz",
        "six",
        "PyQt5",
    ],
)
