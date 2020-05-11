__version__ = '0.5.1'

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
                 connectivity=None,          # Connectivity between atoms
                 connect_thresh=1.10,        # threshold to use in connectivity='auto'
                 permutation=None,           # fix permutation
                 center=None):               # Center of symmetry measure (if None: search optimum)

        conv = _get_connectivity_vector(connectivity, len(coordinates), central_atom)
        if central_atom is None:
            central_atom = 0

        operation = group[0].lower().encode('ascii')
        try:
            operation_axis = int(group[1:])
        except IndexError:
            operation_axis = 1
        except ValueError:
            raise Exception('Wrong symmetry label')

        if center is None:
            fixcenter = False
            center= [0, 0, 0]
        else:
            fixcenter = True

        if connectivity is None:       # not use connectivity
            connect_type = 0
        elif connectivity == 'auto':   # use connectivity from atomic radii
            connect_type = 1
        else:
            connect_type = 2          # custom connectivity by user

        if labels is None:
            labels = ['{}'.format(i) for i in range(len(coordinates))]
            connect_type = 0

        if permutation is None:
            permutation = [0] * len(coordinates)
            fix_permutation = False
        else:
            fix_permutation = True

        labels = np.array([list('{:<2}'.format(char)) for char in labels], dtype='S')
        coordinates = np.ascontiguousarray(coordinates)

        outputs = symgrouplib.symgroup(coordinates, multi, labels, central_atom, operation,
                                       operation_axis, fixcenter, center, connect_type, conv,
                                       connect_thresh, fix_permutation, permutation)

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

    print('---')

    cart_coordinates = [[ 2.756888080, -0.000000000, 0.000000000],
                        [ 2.090404324,  1.404011560, 0.000000000],
                        [ 1.021653372,  2.471437166, 0.000000000],
                        [-0.528304382,  2.333764742, 0.000000000],
                        [-1.813582469,  1.860517986, 0.000000000],
                        [-2.452704063,  0.710363638, 0.000000000],
                        [-2.452704063, -0.710363638, 0.000000000],
                        [-1.813582469, -1.860517986, 0.000000000],
                        [-0.528304382, -2.333764742, 0.000000000],
                        [ 1.021653372, -2.471437166, 0.000000000],
                        [ 2.090404324, -1.404011560, 0.000000000]]

    measure = Symgroupy(coordinates=cart_coordinates,
                        group='c11',
                        labels=['O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O', 'O'],
                        connectivity='auto',
                        )

    print('CSM: {}'.format(measure.csm))
    print('Optimum axis: {}'.format(measure.optimum_axis))
    print('Optimum permutation: {}'.format(measure.optimum_permutation))
    print('Nearest structure')
    print(measure.nearest_structure)
    print('Reference axis')
    print(measure.reference_axis)
    print('multi CMS')
    print(measure.cms_multi)
    print('multi axis')
    print(measure.axis_multi)

    cart_coordinates = [[ 0.506643354, -1.227657970, 0.000000000],
                        [ 1.303068499,  0.000000000, 0.000000000],
                        [ 0.506643354,  1.227657970, 0.000000000],
                        [-0.926250976,  0.939345948, 0.000000000],
                        [-0.926250976, -0.939345948, 0.000000000],
                        ]

    measure = Symgroupy(coordinates=cart_coordinates,
                        group='C5',
                        labels=['C', 'C', 'C', 'C', 'C'],
                        connectivity=[(1,2), (2,3), (3,4), (4,5), (5,1)],
                        )

    print('CSM: {}'.format(measure.csm))
    print('Optimum permutation: {}'.format(measure.optimum_permutation))

    print('*******************')

    measure = Symgroupy(coordinates=cart_coordinates,
                        group='C5',
                        labels=['C', 'C', 'C', 'C', 'C'],
                        connectivity=None,
                        permutation=[5, 1, 2, 3, 4]
                        )

    print('-----------')
    print('CSM: {}'.format(measure.csm))
    print('Optimum axis: {}'.format(measure.optimum_axis))
    print('Optimum permutation: {}'.format(measure.optimum_permutation))
    print('Nearest structure')
    print(measure.nearest_structure)
    print('Reference axis')
    print(measure.reference_axis)
    print('multi CMS')
    print(measure.cms_multi)
    print('multi axis')
    print(measure.axis_multi)
