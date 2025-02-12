***************
Getting Started
***************

Requirements
============

AstroVoigtFit has the following strict dependancies:

* `Python <https://www.python.org/>`_

* `Numpy <https://numpy.org/>`_

* `Astropy <https://www.astropy.org>`_ 4.0 or later

* `Matplotlib <https://matplotlib.org/>`_

* `Pandas <https://pandas.pydata.org/>`_

* `lmfit <https://pypi.org/project/lmfit/>`_

* `specutils <https://pypi.org/project/specutils/>`_

AstroVoigtFit also optionally depends on other packages:

* `Scipy <https://www.scipy.org/>`_ 0.19 or later

* `pytz <https://pypi.org/project/pytz/>`_ .

* `six <https://pypi.org/project/six/>`_ .

* `cycler <https://pypi.org/project/Cycler/>`_ .

* `pyparsing <https://pypi.org/project/pyparsing/>`_ .


For documentation building, AstroVoigtFit depends on `sphinx-astropy
<https://github.com/astropy/sphinx-astropy>`_ (0.4 or later) 

Installation
------------

There are two methods of installing AstroVoigtFit

If you plan to use the AstroVoigtFit package without developing, you can use the regular pip installation::
  
    git clone https://github.com/jancami/AstroVoigtFit.git
    cd edibles
    pip install .

This will install AstroVoigtFit alongside your other Python packages.

If you plan to develop the package, you should use the developer mode of pip installation::

    git clone https://github.com/jancami/AstroVoigtFit.git
    cd edibles
    pip install -e .

AstroVoigtFit should now be setup and you can get to fitting!
