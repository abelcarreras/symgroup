from symgroupy import symgrouplib
import numpy as np


class Symgroupy:
    def __init__(self,
                 coordinates,
                 multi,
                 group,
                 labels=None,
                 central_atom=0):

        oper = b'{}'.format(group[0].lower())
        oeix = int(group[1])
        pm = 1
        nll = 4

        labels = np.array([list('{:<2}'.format(char)) for char in labels], dtype='S')
        coordinates = np.array(coordinates)

        symgrouplib.symgroup('test', coordinates, multi, labels, central_atom, oper, oeix)

        print('test')

if __name__ == '__main__':

    #symgrouplib.symgroup('test')

    coordinates = [[5.51829, -1.68040, 22.81703],
                   [6.78978, -3.22298, 23.08474],
                   [6.27712, -0.12712, 21.76775],
                   [4.24692, -3.22298, 22.54931],
                   [4.75958, -0.12712, 23.86630]]

    fen4 = Symgroupy(coordinates=coordinates,
                     multi=8,
                     group='C3',
                     labels=['Fe', 'N', 'N', 'N', 'N'],
                     central_atom=1)
