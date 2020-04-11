__version__ = '0.4.1'

from symgroupy import symgrouplib
import numpy as np


def _get_connectivity_vector(connect, na, central_atom=0):

    matrix = np.zeros((na, na), dtype=float)

    try:
        for pair in connect:
            matrix[pair[0]-1, pair[1]-1] = 1.0
            matrix[pair[1]-1, pair[0]-1] = 1.0

        if central_atom > 0:
            matrix = np.vstack((matrix, matrix[central_atom-1, :]))
            matrix = np.hstack((matrix, matrix[:, central_atom-1][None].T))

            matrix = np.delete(matrix, central_atom-1, axis=0)
            matrix = np.delete(matrix, central_atom-1, axis=1)
    except TypeError:
        pass

    return np.ascontiguousarray(matrix[np.tril_indices(na)])


class Symgroupy:
    def __init__(self,
                 coordinates,                # coordinates in angstrom
                 group,                      # symmetry point group
                 multi=1,                    # multiple measures
                 labels=None,                # atomic symbols (or other representative labels)
                 central_atom=None,          # Atom number that contains the center atom (if exist)
                 connectivity='auto',        # use connectivity from atomic radii
                 connect_thresh=1.10,        # threshold to use in connectivity='auto'
                 fix_permutation=False,      # fix permutation
                 center=None):               # Center of symmetry measure (if None: search optimum)

        conv = _get_connectivity_vector(connectivity, len(coordinates), central_atom)
        if central_atom is None:
            central_atom = 0

        operation = group[0].lower().encode('ascii')
        try:
            operation_axis = int(group[1])
        except IndexError:
            operation_axis = 1

        if center is None:
            fixcenter = False
            center= [0, 0, 0]
        else:
            fixcenter = True

        if connectivity is None:
            connect_type = 0
        elif connectivity == 'auto':
            connect_type = 1
        else:
            connect_type = 2

        if labels is None:
            labels = ['{}'.format(i) for i in range(len(coordinates))]
            connect_type = 0

        labels = np.array([list('{:<2}'.format(char)) for char in labels], dtype='S')
        coordinates = np.ascontiguousarray(coordinates)

        outputs = symgrouplib.symgroup(coordinates, multi, labels, central_atom, operation,
                                       operation_axis, fixcenter, center, connect_type, conv,
                                       connect_thresh, fix_permutation)

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

    cart_coordinates = [[-0.15936255,   -0.27888446,    0.00000000],
                        [0.19729187,   -1.28769446,    0.00000000],
                        [0.19731029,    0.22551373,    0.87365150],
                        [0.19731029,    0.22551373,   -0.87365150],
                        [-1.22936255,   -0.27887127,    0.00000000]]


# Symmetry point groups available:
    # e: identity
    # i: inversion
    # r: reflection
    # cn: rotation (n:order)
    # sn: improper rotation (n:order)

    bonds = [(1, 2), (1, 3), (1, 4), (1, 5)]  # pairs of atoms

    fen4 = Symgroupy(coordinates=cart_coordinates,
                     group='c4',
                     multi=3,
                     labels=['C', 'H', 'H', 'H', 'H'],
                     central_atom=1,
                     fix_permutation=False,
                     center=[-0.15936255,   -0.27888446,    0.00000000],
                     connectivity='auto',
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
