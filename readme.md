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
   ```shell
   ./configure (see --help for available options)
   cd build
   make install
   ```
2b. Compile as a python module
   ```shell
   cd python
   python setup.py install --user
   ```
2c. Obtain from PyPi repository
```shell
pip install symgroupy --user
```
Python API
----------
```python
from symgroupy import Symgroupy

fen4 = Symgroupy(coordinates=[[15.5182, -1.68040, 22.81703],
                              [6.78978, -3.22298, 23.08474],
                              [6.27712, -0.12712, 21.76775],
                              [4.24692, -3.22298, 22.54931],
                              [4.75958, -0.12712, 23.86630]],
                 group='c3',
                 multi=8,
                 labels=['Fe', 'N', 'N', 'N', 'N'],
                 central_atom=1)

print('CSM: {}'.format(fen4.csm))
print('Optimum axis: {}'.format(fen4.optimum_axis))
print('Optimum permutation: {}'.format(fen4.optimum_permutation))
print('Nearest structure')
print(fen4.nearest_structure)
print('Reference axis')
print(fen4.reference_axis)
print('multi CMS')
print(fen4.cms_multi)
print('multi axis')
print(fen4.axis_multi)

```

Authors
--------------------------------------------------------

This software has been developed by David Casanova
<br>Python module by Abel Carreras

The theoretical background implemented in this software is described in:
<br>Pinsky M, Dryzun C, Casanova D, Alemany P, Avnir D, J Comput Chem. 29:2712-21 (2008)
<br>Pinsky M, Casanova D, Alemany P, Alvarez S, Avnir D, Dryzun C, Kizner Z, Sterkin A.  J Comput Chem. 29:190-7 (2008)