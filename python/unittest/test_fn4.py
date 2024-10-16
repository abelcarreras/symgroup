from symgroupy import Symgroupy
import numpy as np
import unittest


class TestFeN4(unittest.TestCase):

    def setUp(self):

        cart_coordinates = [[5.51829, -1.68040, 22.81703],
                            [6.78978, -3.22298, 23.08474],
                            [6.27712, -0.12712, 21.76775],
                            [4.24692, -3.22298, 22.54931],
                            [4.75958, -0.12712, 23.86630]]

        self.fen4 = Symgroupy(coordinates=cart_coordinates,
                              group='c3',
                              multi=8,
                              labels=['Fe', 'N', 'N', 'N', 'N'],
                              central_atom=1,
                              connectivity=None,
                              )

    def test_csm(self):
        csm_ref = 8.70279043841
        print(self.fen4.csm)
        np.testing.assert_almost_equal(self.fen4.csm, csm_ref, decimal=8)

    def test_optimum_axis(self):
            optimum_axis_ref = [0.74666007, 0.57751368, 0.33011617]

            np.testing.assert_allclose(optimum_axis_ref, self.fen4.optimum_axis, rtol=1e-6)

    def test_optimum_permutation(self):
        permutation_ref = [1, 5, 2, 4, 3]
        np.testing.assert_allclose(permutation_ref, self.fen4.optimum_permutation, rtol=1e-6)

    def test_nearest_structure(self):
        nearest_structure_test = [[ 5.518322,   -1.67754667, 22.81702733],
                                  [ 6.9599505,  -2.78838476, 23.4545737 ],
                                  [ 6.15622621, -0.56178074, 21.37716726],
                                  [ 4.07651706, -2.79131466, 22.17956268],
                                  [ 4.88252957, -0.56155244, 24.25762773]]

        np.testing.assert_allclose(nearest_structure_test, self.fen4.nearest_structure, rtol=1e-3, atol=1e-2)

    def test_reference_axis(self):

        reference_axis_ref = [[ 0.00000000,  0.00000000,  0.00000000],
                              [ 0.23511513,  0.23511513, -0.94310218],
                              [-0.62226969,  0.78179204,  0.03976887],
                              [ 0.74666007,  0.57751368,  0.33011617]]

        np.testing.assert_allclose(reference_axis_ref, self.fen4.reference_axis, rtol=1e-6)

    def test_multi_csm(self):

        multi_csm_8_ref = [ 8.70279044,  8.70279106,  8.70279249,  8.70279305,
                            41.25136234, 63.29000792, 63.29001106, 63.42240444]

        np.testing.assert_allclose(multi_csm_8_ref, self.fen4.csm_multi, rtol=1e-6)

    def test_multi_axis(self):

        multi_axis_8_ref = [[ 7.46660070e-01,  5.77513683e-01,  3.30116172e-01],
                            [ 3.30209500e-01,  5.77187379e-01, -7.46871084e-01],
                            [-3.30209474e-01,  5.77187371e-01,  7.46871102e-01],
                            [ 7.46660071e-01, -5.77513670e-01,  3.30116192e-01],
                            [ 5.10308545e-08,  1.00000000e+00, -9.38206158e-09],
                            [ 6.85623259e-01,  7.27918284e-01,  7.46450811e-03],
                            [ 6.85623255e-01, -7.27918287e-01,  7.46451763e-03],
                            [ 5.02329182e-01,  7.29481935e-01, -4.64242931e-01]]

        np.testing.assert_allclose(multi_axis_8_ref, self.fen4.axis_multi, rtol=1e-6)



class TestFeN4_2(unittest.TestCase):

    def setUp(self):

        cart_coordinates = [[ 6.78978, -3.22298, 23.08474],
                            [ 6.27712, -0.12712, 21.76775],
                            [15.51829, -1.68040, 22.81703],
                            [ 4.24692, -3.22298, 22.54931],
                            [ 4.75958, -0.12712, 23.86630]]

        self.fen4 = Symgroupy(coordinates=cart_coordinates,
                              group='c3',
                              multi=3,
                              labels=['N', 'N', 'Fe', 'N', 'N'],
                              central_atom=3,
                              center=[0, 0, 0],
                              connectivity='auto',
                              )

    def test_csm(self):
        csm_ref = 2.631096544254
        print(self.fen4.csm)
        np.testing.assert_almost_equal(self.fen4.csm, csm_ref, decimal=8)

    def test_optimum_axis(self):
        optimum_axis_ref = [0.3211507,  -0.07036835,  0.94441015]

        np.testing.assert_allclose(optimum_axis_ref, self.fen4.optimum_axis, rtol=1e-6)

    def test_optimum_permutation(self):
        permutation_ref = [1, 4, 3, 5, 2]
        np.testing.assert_allclose(permutation_ref, self.fen4.optimum_permutation, rtol=1e-6)

    def test_reference_axis(self):

        reference_axis_ref = [[ 0.00000000,  0.00000000,  0.00000000],
                              [ 0.69496185,  0.69496185, -0.18454288],
                              [-0.64334303,  0.71559507,  0.27209079],
                              [ 0.32115070, -0.07036835,  0.94441015]]

        np.testing.assert_allclose(reference_axis_ref, self.fen4.reference_axis, rtol=1e-6)

    def test_multi_csm(self):

        multi_csm_8_ref = [ 2.631097, 2.879698, 2.890738]

        np.testing.assert_allclose(multi_csm_8_ref, self.fen4.csm_multi, rtol=1e-6)

    def test_multi_axis(self):

        multi_axis_3_ref = [[0.32115070, -0.07036835,  0.94441015],
                            [0.32009532, -0.06959149,  0.94482592],
                            [0.32026730, -0.06985121,  0.94474847]]

        np.testing.assert_allclose(multi_axis_3_ref, self.fen4.axis_multi, rtol=1e-6)


class TestFeN4_groups(unittest.TestCase):

    def setUp(self):

        self._cart_coordinates = [[5.51829, -1.68040, 22.81703],
                                  [6.78978, -3.22298, 23.08474],
                                  [6.27712, -0.12712, 21.76775],
                                  [4.24692, -3.22298, 22.54931],
                                  [4.75958, -0.12712, 23.86630]]

    def test_cs(self):

        fen4 = Symgroupy(coordinates=self._cart_coordinates,
                         group='cs',
                         labels=['Fe', 'N', 'N', 'N', 'N'],
                         central_atom=1,
                         connectivity=None,
                         )

        csm_ref = 1.7804508
        print(fen4.csm)
        np.testing.assert_almost_equal(fen4.csm, csm_ref, decimal=3)

    def test_ci(self):

        fen4 = Symgroupy(coordinates=self._cart_coordinates,
                         group='ci',
                         labels=['Fe', 'N', 'N', 'N', 'N'],
                         central_atom=1,
                         connectivity=None,
                         )

        csm_ref = 12.24274
        print(fen4.csm)
        np.testing.assert_almost_equal(fen4.csm, csm_ref, decimal=3)

    def test_s1(self):

        fen4 = Symgroupy(coordinates=self._cart_coordinates,
                         group='cs',
                         multi=8,
                         labels=['Fe', 'N', 'N', 'N', 'N'],
                         central_atom=1,
                         connectivity=None,
                         )

        csm_ref = 1.78045
        print(fen4.csm)
        np.testing.assert_almost_equal(fen4.csm, csm_ref, decimal=3)

    def test_s2(self):

        fen4 = Symgroupy(coordinates=self._cart_coordinates,
                         group='s2',
                         multi=8,
                         labels=['Fe', 'N', 'N', 'N', 'N'],
                         central_atom=1,
                         connectivity=None,
                         )

        csm_ref = 12.24274
        print(fen4.csm)
        np.testing.assert_almost_equal(fen4.csm, csm_ref, decimal=3)
