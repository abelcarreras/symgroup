Example
=======

Simple input data file (**.zdat**) ::

 $Ge9 rotation analysis
 9 0
 c 4
 BEHLOQ
 Ge2    25.05753   17.90064   15.10444
 Ge5    24.55283   20.39220   13.65426
 Ge3    22.93143   16.75130   13.88961
 Ge1    22.77389   19.09191   15.02517
 Ge8    25.39808   16.13914   13.16510
 Ge9    26.22494   18.55152   12.87217
 Ge4    22.43074   19.19254   12.44356
 Ge6    24.67853   19.59012   11.11344
 Ge7    23.82786   17.18736   11.46072
 BEHLUW
 Ge2     6.07390    8.45936   14.02226
 Ge5     8.83888    8.96150   14.18831
 Ge3     6.15319    9.10288   11.33143
 Ge1     7.14830   10.62999   13.14511
 Ge8     5.85402    6.63872   12.18891
 Ge9     7.86582    6.52045   13.78435
 Ge4     8.88800    9.59116   11.54514
 Ge6     9.90916    7.34901   12.47335
 Ge7     7.92895    7.26098   10.83291

Output file (**.ztab**) ::

 Ge9 rotation analysis
 3D Rotation Measures C 4 Group

 STRUCTURE         MEASURE
  BEHLOQ     ,     0.58207
  BEHLUW     ,     0.69466

Input file with keywords::

 $ ML4 Structures with d6 configuration %label
 %multi 8
 4 1
 c 3
 CIWQUU
 Fe    5.51829   -1.68040   22.81703
 N     6.78978   -3.22298   23.08474
 N     6.27712   -0.12712   21.76775
 N     4.24692   -3.22298   22.54931
 N     4.75958   -0.12712   23.86630

In this example, the program searches for the best eight **C:sub:`3`** axis and
compares only atoms with the same label.

Output file (**.zout2**)::

 ML4 Structures with d6 configuration

 ************************************************
 CIWQUU

  6.52709   6.52709   6.52709   6.52709   6.52709

 -0.74666   0.74666  -0.33021   0.33021  -0.33021
 -0.57751   0.57751  -0.57719   0.57719   0.57719
 -0.33012   0.33012   0.74687  -0.74687   0.74687

  6.52709   6.52709   6.52709

  0.33021  -0.74666   0.74666
 -0.57719   0.57751  -0.57751
 -0.74687  -0.33012   0.33012

The **.zout2** file contains the 8 lowest :math:`a^2 + b^2 = c^2` measures with the corresponding axis orientations.