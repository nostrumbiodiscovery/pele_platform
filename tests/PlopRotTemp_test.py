import unittest
from PlopRotTemp import intersect_tors, add_tors, remove_tors, find_names_in_mae

class TestPlopRotTemp(unittest.TestCase):

    def test_intersectTors(self):
        torsion1 = [[1, 3], [7, 9], [10, 15], [11, 12]]
        torsion2 = [[1, 2], [7, 9], [11, 12]]
        correct_intersection = [[7, 9], [11, 12]]
        intersect = intersect_tors(torsion1, torsion2)
        self.assertEqual(intersect, correct_intersection)

    def test_addtors(self):
        torsion1 = [[1, 3], [7, 9], [10, 15], [11, 12]]
        torsion2 = [[1, 2], [7, 9], [11, 12]]
        correct_addition = [[1, 3], [7, 9], [10, 15], [11, 12], [1, 2]]
        addition = add_tors(torsion1, torsion2)
        self.assertEqual(addition, correct_addition)

    def test_removetors(self):
        torsion1 = [[1, 3], [7, 9], [10, 15], [11, 12]]
        torsion2 = [[1, 2], [7, 9], [11, 12]]
        correct_substraction1_2 = [[1, 3], [10, 15]]
        subs1_2 = remove_tors(torsion1, torsion2)
        correct_substraction2_1 = [[1, 2]]
        subs2_1 = remove_tors(torsion2, torsion1)
        self.assertEqual(subs1_2, correct_substraction1_2)
        self.assertEqual(subs2_1, correct_substraction2_1)

    def test_find_names_in_mae(self):
        filename_2012v = "tests_data/2gp_2012v.mae"
        filename_2016v = "tests_data/2gp_2016v.mae"
        atom_names = [' P1 ', ' O1 ', ' O2 ', ' O3 ', ' O4 ', ' C1 ', ' C2 ', ' O5 ', ' C3 ', ' O6 ', ' C4 ', ' O7 ',
                      ' C5 ', ' N1 ', ' C6 ', ' N2 ', ' C7 ', ' C8 ', ' O8 ', ' N3 ', ' C9 ', ' N4 ', ' N5 ', ' C10',
                      ' H1 ', ' H2 ', ' H3 ', ' H4 ', ' H5 ', ' H6 ', ' H7 ', ' H8 ', ' H9 ', ' H10', ' H11', ' H12']
        names_2012 = find_names_in_mae(filename_2012v)
        names_2016 = find_names_in_mae(filename_2016v)
        self.assertEqual(names_2012, atom_names)
        self.assertEqual(names_2016, atom_names)
