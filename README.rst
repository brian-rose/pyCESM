================
pyCESM
================
----------
 Python package for interactive work with the Community Earth System Model
----------

Author
--------------
| **Brian E. J. Rose**
| Department of Atmospheric and Environmental Sciences
| University at Albany
| brose@albany.edu

Installation
----------------
``python setup.py``

    or, if you are developing new code

``python setup.py develop``


About pyCESM
--------------
Currently consists of two parts:

- Python / numpy implementations of some of the 
  thermodynamic routines (e.g. saturation vapor pressure and specific humidity)
  used by the CAM atmospheric model.

- A Python interface for working with CAM output that automates some grid-aware diagnostic calculations. 
  This functionality is built on top of xarray_.

Based on (and tested with) version 1.2.1 of CESM

License
---------------
This code is freely available under the MIT license.
See the accompanying LICENSE file.

.. _xarray: http://xarray.pydata.org
