How to install
==============

Standalone binary
-----------------
Use configure script located in the symgroup directory to create
a suitable makefile with the available compiler options.

   $ ./configure [options]

Options:
  release: Use strong optimizations (slow to compile fast to run)
  debug: Use few code optimizations (fast to compile slow to run)
  gcc : Use gfortran compiler (default)
  intel: Use Intel fortran compiler (ifort)
  mkl: Use Intel Math kernel Library ::

   $ cd build
   $ make symgroup

the executable will be placed in EXE directory

Fortran library
---------------
To build symgroup as a fortran library use configure script as described above
and type::

   $ cd build
   $ make symgrouplib

the library will be placed in LIB directory

Python module (symgroupy)
-------------------------

Python module for symgroup is called symgroupy and it is compatible with
both python2 and python3. To install symgroupy go to python directory
and use :file: `setup.py` to install the code using :program: `setuptools` python
module. A simple setup may be ::

   $ setup.py install --user

the code will be installed as a python module. To check that it is properly installed you can
run the :program: `python interpret and execute ::

   import symgroupy

if the module should be loaded without warnings and errors