
Input file
==========

STRUCTURE
---------
1. One line with user information to identify the case under study (optional). If this line is
present, it must have a $ symbol in the first column.
2. Keywords for the options chosen (optional).
3. This line must contain two integer numbers that indicate the number of atoms (disregarding the central atom) and the position of the central atom in the list of atomic coordinates. If a central atom is absent, the second number must be a 0.
4. Symmetry element to be measured: E for the identity, i for inversion, r for a reflection,
c n for a proper rotation, where n is the order of the rotation (i.e., c 2, or c 3, or...), and s n for an improper rotation.
5. Label identifying the compound to which the subsequent list of atomic coordinates belongs. This might be the CSD refcode, the ICSD collection code, or the compound name.
6. Atomic label and coordinates, one atom per line (free format). If the %connect keyword is used (see below), the label must correspond to an element symbol. If no atomic labels are given, the %nolabel keyword must be used (see below).
A set of lines 5 and 6 can be added for as many atom sets as desired.

NOTES
-----
- Type SYMGROUP + to obtain a list of the symmetry elements that the program can measure.
- Type SYMGROUP +keywords to obtain a list of the available keywords on the screen.
- Blank lines can be introduced anywhere in the input file. Comment lines, indicated by the ! symbol in the first column can be added at any point of the input file.
- The program is case insensitive, upper or lower case letters can be used for symmetry element labels or keywords in the .zdat file.