"""
module to test data_types
"""

import unittest

import pathadd
from silica_lattice_model import data_types as dt


class TestBccLattice3DN1list(unittest.TestCase):
    """ class to test BccLattice3D """

    def test_single_unit(self):
        lat_obj = dt.BccLattice3D(2, 2, 2)
        # lat_obj = dt.BccLattice3D(1, 1, 1)
        self.assertEqual(lat_obj.nsites, 16)

        # test_data_i = (ni, irx, iry, irz)
        test_data = (
            (0, 2, 2, 2),
            (1, 4, 2, 2),
            (2, 4, 4, 2),
            (3, 2, 4, 2),
            (4, 2, 2, 4),
            (5, 4, 2, 4),
            (6, 4, 4, 4),
            (7, 2, 4, 4),
        )

        for isite in range(lat_obj.nsites):
            isite_rx = lat_obj.rx[isite]
            isite_ry = lat_obj.ry[isite]
            isite_rz = lat_obj.rz[isite]

            if isite_rx == 3 and isite_ry == 3 and isite_rz == 3:

                for idx, correct_rx, correct_ry, correct_rz in test_data:
                    irx = lat_obj.rx[lat_obj.n1list[isite, idx]]
                    iry = lat_obj.ry[lat_obj.n1list[isite, idx]]
                    irz = lat_obj.rz[lat_obj.n1list[isite, idx]]

                    # print('irx', irx, iry, irz)
                    # print('correct_irx', correct_rx, correct_ry, correct_rz)
                    # print()

                    self.assertEqual(irx, correct_rx)
                    self.assertEqual(iry, correct_ry)
                    self.assertEqual(irz, correct_rz)
        pass


class TestBccLattice3DN2list(unittest.TestCase):
    """ class to test BccLattice3D """

    def test_single_unit(self):
        lat_obj = dt.BccLattice3D(5, 5, 5)
        # lat_obj = dt.BccLattice3D(1, 1, 1)

        test_data = (
            (0, 1, 3, 3),
            (1, 5, 3, 3),
            (2, 3, 1, 3),
            (3, 3, 5, 3),
            (4, 3, 3, 1),
            (5, 3, 3, 5)
        )
        for isite in range(lat_obj.nsites):
            isite_rx = lat_obj.rx[isite]
            isite_ry = lat_obj.ry[isite]
            isite_rz = lat_obj.rz[isite]

            if isite_rx == 3 and isite_ry == 3 and isite_rz == 3:

                for idx, correct_rx, correct_ry, correct_rz in test_data:
                    irx = lat_obj.rx[lat_obj.n2list[isite, idx]]
                    iry = lat_obj.ry[lat_obj.n2list[isite, idx]]
                    irz = lat_obj.rz[lat_obj.n2list[isite, idx]]

                    # print('irx', irx, iry, irz)
                    # print('correct_irx', correct_rx, correct_ry, correct_rz)
                    # print()

                    self.assertEqual(irx, correct_rx)
                    self.assertEqual(iry, correct_ry)
                    self.assertEqual(irz, correct_rz)
        pass
