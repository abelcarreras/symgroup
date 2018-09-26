[![Build Status](https://travis-ci.org/abelcarreras/symgroup.svg?branch=master)](https://travis-ci.org/abelcarreras/symgroup)
[![Coverage Status](https://coveralls.io/repos/github/abelcarreras/symgroup/badge.svg?branch=master)](https://coveralls.io/github/abelcarreras/symgroup?branch=master)
[![PyPI version](https://badge.fury.io/py/symgroupy.svg)](https://pypi.org/project/symgroupy)

SymGroup
=========
Software to calculate continuous symmetry measures of 
molecular structures


Installation instructions
---------------------------------------------------------

1. Requirements
  - Lapack & Blas libraries
  - Fortran77 compiler
  - cmake 2.6
  - (optional) MKL


2a. Install as standalone binary
   ```
   ./configure (see --help for available options)
   cd build
   make install
   ```
2b. Compile as a python module
   ```
   cd python
   python setup.py install --user
   ```
2c. Obtain from PyPi repository
```
pip install symgroupy --user
```

Authors
--------------------------------------------------------

This software has been developed by David Casanova
<br>Python module by Abel Carreras