from symgroupy import symgrouplib
import numpy as np


class Symgroupy:
    def __init__(self,
                 coordinates,
                 group,
                 multi=1,
                 labels=None,
                 central_atom=None):

        if central_atom is None:
            central_atom_index = 0
        else:
            central_atom_index = central_atom

        operation = b'{}'.format(group[0].lower())
        try:
            operation_axis = int(group[1])
        except IndexError:
            operation_axis = 1

        labels = np.array([list('{:<2}'.format(char)) for char in labels], dtype='S')
        coordinates = np.array(coordinates, dtype='double', order='c')

        outputs = symgrouplib.symgroup(coordinates, multi, labels, central_atom_index, operation, operation_axis)

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

    #symgrouplib.symgroup('test')

    coordinates = [[5.51829, -1.68040, 22.81703],
                   [6.78978, -3.22298, 23.08474],
                   [6.27712, -0.12712, 21.76775],
                   [4.24692, -3.22298, 22.54931],
                   [4.75958, -0.12712, 23.86630]]

    fen4 = Symgroupy(coordinates=coordinates,
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