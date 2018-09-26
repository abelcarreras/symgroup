EXAMPLE
=======

Simple input data file (.zdat) ::

 $Ge9 rotation analysis
 9 0
 c4
  BEHLOQ Ge2 Ge5 Ge3 Ge1 Ge8 Ge9 Ge4 Ge6 Ge7 BEHLUW Ge2 Ge5 Ge3 Ge1 Ge8 Ge9 Ge4 Ge6 Ge7
Output file (.ztab) Ge9 rotation
3D Rotation
STRUCTURE
BEHLOQ , 0.58207 BEHLUW , 0.69466
25.05753 17.90064 24.55283 20.39220 22.93143 16.75130 22.77389 19.09191 25.39808 16.13914 26.22494 18.55152 22.43074 19.19254 24.67853 19.59012 23.82786 17.18736
6.07390 8.45936 8.83888 8.96150 6.15319 9.10288 7.14830 10.62999 5.85402 6.63872 7.86582 6.52045 8.88800 9.59116 9.90916 7.34901 7.92895 7.26098
analysis
Measures C 4 Group
15.10444
13.65426
13.88961
15.02517
13.16510
12.87217
12.44356
11.11344
11.46072
14.02226
14.18831
11.33143
13.14511
12.18891
13.78435
11.54514
12.47335
10.83291
 MEASURE
Input file with keywords:
$ ML4 Structures with d6 configuration %label
%multi 8
41
c3
CIWQUU
Fe 5.51829 -1.68040 22.81703 N 6.78978 -3.22298 23.08474 N 6.27712 -0.12712 21.76775 N 4.24692 -3.22298 22.54931 N 4.75958 -0.12712 23.86630
In this example, the program searches for the best eight C3 axis and compares only atoms with the same label.
Output file (.zout2)
ML4 Structures with d6 configuration
************************************** CIWQUU
  6.52709 6.52709 6.52709
-0.74666 0.74666 -0.33021 -0.57751 0.57751 -0.57719 -0.33012 0.33012 0.74687
6.52709 6.52709 6.52709
0.33021 -0.74666 0.74666 -0.57719 0.57751 -0.57751 -0.74687 -0.33012 0.33012
6.52709 6.52709
0.33021 -0.33021
0.57719 0.57719 -0.74687 0.74687
**************************************
The . zout2 file contains the 8 lowest C3 measures with the corresponding axis orientations.