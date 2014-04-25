#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for Screw6 class."""


import unittest

from sympy import Matrix
from sympy import ShapeError
from sympy import zeros

from pysymoro.screw6 import Screw6


class TestScrew6(unittest.TestCase):
    """Unit test for Screw6 class."""
    def setUp(self):
        # setup data matrices to compare against
        self.data = Matrix([
            [1, 2, 3, 4, 5, 6],
            [7, 8, 9, 10, 11, 12],
            [13, 14, 15, 16, 17, 18],
            [19, 20, 21, 22, 23, 24],
            [25, 26, 27, 28, 29, 30],
            [31, 32, 33, 34, 35, 36]
        ])
        self.data_tl = Matrix([
            [1, 2, 3],
            [7, 8, 9],
            [13, 14, 15]
        ])
        self.data_tr = Matrix([
            [4, 5, 6],
            [10, 11, 12],
            [16, 17, 18]
        ])
        self.data_bl = Matrix([
            [19, 20, 21],
            [25, 26, 27],
            [31, 32, 33]
        ])
        self.data_br = Matrix([
            [22, 23, 24],
            [28, 29, 30],
            [34, 35, 36]
        ])
        # setup Screw6 instances
        self.empty = Screw6()
        self.indiv = Screw6()
        self.indiv.val = self.data

    def test_init(self):
        """Test constructor."""
        # test instance type
        self.assertIsInstance(self.empty, Screw6)
        self.assertIsInstance(self.indiv, Screw6)
        self.assertIsInstance(Screw6(self.data), Screw6)
        self.assertIsInstance(Screw6(value=self.data), Screw6)
        self.assertIsInstance(
            Screw6(
                self.data_tl, self.data_tr,
                self.data_bl, self.data_br
            ), Screw6
        )
        self.assertIsInstance(
            Screw6(
                tl=self.data_tl, bl=self.data_bl,
                tr=self.data_tr, br=self.data_br
            ), Screw6
        )
        # test raise appropriate exception
        self.assertRaises(NotImplementedError, Screw6, 3, 4)
        self.assertRaises(NotImplementedError, Screw6, tl=5, tr=6)

    def test_val(self):
        """Test get and set of val()"""
        # test get
        self.assertEqual(self.empty.val, zeros(6, 6))
        self.assertEqual(self.indiv.val, self.data)
        # test set
        self.indiv.val = self.data.transpose()
        self.assertEqual(self.indiv.val, self.data.transpose())
        with self.assertRaises(ShapeError):
            self.empty.val = Matrix([3, 3])

    def test_topleft(self):
        """Test get and set of topleft()"""
        # test get
        self.assertEqual(self.empty.topleft, zeros(3, 3))
        # test set
        self.indiv.topleft = self.data_tl
        self.assertEqual(self.indiv.topleft, self.data_tl)
        with self.assertRaises(ShapeError):
            self.empty.topleft = Matrix([3, 3])

    def test_topright(self):
        """Test get and set of topright()"""
        # test get
        self.assertEqual(self.empty.topright, zeros(3, 3))
        # test set
        self.indiv.topright = self.data_tr
        self.assertEqual(self.indiv.topright, self.data_tr)
        with self.assertRaises(ShapeError):
            self.empty.topright = Matrix([3, 3])

    def test_botleft(self):
        """Test get and set of botleft()"""
        # test get
        self.assertEqual(self.empty.botleft, zeros(3, 3))
        # test set
        self.indiv.botleft = self.data_bl
        self.assertEqual(self.indiv.botleft, self.data_bl)
        with self.assertRaises(ShapeError):
            self.empty.botleft = Matrix([3, 3])

    def test_botright(self):
        """Test get and set of botright()"""
        # test get
        self.assertEqual(self.empty.botright, zeros(3, 3))
        # test set
        self.indiv.botright = self.data_br
        self.assertEqual(self.indiv.botright, self.data_br)
        with self.assertRaises(ShapeError):
            self.empty.botright = Matrix([3, 3])

    def test_equality(self):
        """Test __eq__() and __ne__()"""
        self.assertEqual(self.empty, Screw6())
        self.assertNotEqual(self.indiv, Screw6())
        self.assertRaises(
            ValueError, self.empty.__eq__,
            Matrix([1, 2, 3, 4, 5, 6, 7, 8])
        )
        self.assertRaises(
            ValueError, self.empty.__ne__,
            Matrix([1, 2, 3, 4, 5, 6, 7, 8])
        )

def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(TestScrew6)
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


