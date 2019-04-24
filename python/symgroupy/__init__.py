__version__ = '0.2.2'

from symgroupy import symgrouplib
import numpy as np


class Symgroupy:
    def __init__(self,
                 coordinates,  # Cartesian coordinates
                 group,        # symmetry point group
                 multi=1,      # multiple measures
                 labels=None,  # atomic symbols (or other representative labels)
                 central_atom=None):  # Atom number that contains the center atom (if exist)

        if central_atom is None:
            central_atom = 0

        operation = group[0].lower().encode('ascii')
        try:
            operation_axis = int(group[1])
        except IndexError:
            operation_axis = 1

        if labels is None:
            labels = ['{}'.format(i) for i in range(len(coordinates))]

        labels = np.array([list('{:<2}'.format(char)) for char in labels], dtype='S')
        coordinates = np.ascontiguousarray(coordinates)

        outputs = symgrouplib.symgroup(coordinates, multi, labels, central_atom, operation, operation_axis)

        # Reorganize outputs
        self._csm = outputs[0]
        self._nearest_structure = outputs[1]
        self._optimum_axis = outputs[2]
        self._optimum_permutation = outputs[3]
        self._reference_axis = outputs[4]
        self._csm_multi = outputs[5][:multi]
        self._axis_multi = outputs[6][:multi,:]

    @property
    def csm(self):
        return self._csm

    @property
    def nearest_structure(self):
        return self._nearest_structure

    @property
    def optimum_axis(self):
        return self._optimum_axis

    @property
    def optimum_permutation(self):
        return self._optimum_permutation

    @property
    def reference_axis(self):
        return self._reference_axis

    @property
    def cms_multi(self):
        return self._csm_multi

    @property
    def axis_multi(self):
        return self._axis_multi


if __name__ == '__main__':

    cart_coordinates = [[1.51829, -1.68040, 22.81703],
                        [6.78978, -3.22298, 23.08474],
                        [6.27712, -0.12712, 21.76775],
                        [4.24692, -3.22298, 22.54931],
                        [4.75958, -0.12712, 23.86630]]

# Symmetry point groups available:
    # e: identity
    # i: inversion
    # r: reflection
    # cn: rotation (n:order)
    # sn: improper rotation (n:order)

    fen4 = Symgroupy(coordinates=cart_coordinates,
                     group='c3',
                     multi=8,
                     labels=['Fe', 'N', 'N', 'N', 'N'],
                     central_atom=1
                     )

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