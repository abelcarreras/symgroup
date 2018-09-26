.. highlight:: rst

========================
Introduction to symgroup
========================

This program calculates the Continuous Symmetry Measures (CSM), of one or more sets of atomic
coordinates relative to a cyclic point group. The measures can adopt values between 0 and 100.
A zero value indicates that the set of coordinates contains the symmetry element.
The input data file must have the extension .zdat and the results are written in two files with
the same name as the input file and extensions .ztab and .zout. The .ztab file contains the value
of the symmetry measures for each structure and the .zout file contains additional information,
such as the transformed coordinates of the set of atoms and the orientation of the axis associated
to the studied symmetry operation.
The program is invoked by the instruction ::

    $./symgroup inputfile

The input data is read from the inputfile.zdat file, while symgroup is the name
with which the executable version of the program is stored.

.. toctree::
   :maxdepth: 2
   :caption: Table of contents

   install
   input
   keywords