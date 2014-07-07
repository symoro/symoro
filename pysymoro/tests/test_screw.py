#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for Screw class."""


import unittest

from sympy import Matrix
from sympy import ShapeError
from sympy import zeros

from pysymoro.screw import Screw


class TestScrew(unittest.TestCase):
    """Unit test for Screw class."""
    def setUp(self):
        # screw with 0s
        self.empty = Screw()
        # screw created by separate linear and angular terms
        self.indiv = Screw(lin=Matrix([1, 2, 3]), ang=Matrix([4, 5, 6]))
        # screw created by 6x1 matrix
        #self.full = Screw(Matrix([6, 5, 4, 3, 2, 1]))

    def test_init(self):
        """Test constructors."""
        # test instance type
        self.assertIsInstance(self.empty, Screw)
        self.assertIsInstance(self.indiv, Screw)
        # test raise appropriate exception
        self.assertRaises(
            ShapeError, Screw, Matrix([1, 2]), Matrix([7, 8])
        )
        self.assertRaises(
            ShapeError, Screw, Matrix([1, 2, 3, 4, 5, 6, 7, 8])
        )

    def test_val(self):
        """Test get and set of val()"""
        # test get
        self.assertEqual(self.empty.val, zeros(6, 1))
        self.assertEqual(self.indiv.val, Matrix([1, 2, 3, 4, 5, 6]))
        # test set
        self.indiv.val = Matrix([6, 5, 4, 3, 2, 1])
        self.assertEqual(self.indiv.val, Matrix([6, 5, 4, 3, 2, 1]))
        with self.assertRaises(ShapeError):
            self.empty.val = Matrix([3, 3])

    def test_lin(self):
        """Test get and set of lin()"""
        # test get
        self.assertEqual(self.empty.lin, zeros(3, 1))
        # test set
        self.indiv.lin = Matrix([1, 2, 3])
        self.assertEqual(self.indiv.lin, Matrix([1, 2, 3]))
        with self.assertRaises(ShapeError):
            self.empty.lin = Matrix([3, 3])

    def test_ang(self):
        """Test get and set of ang()"""
        # test get
        self.assertEqual(self.empty.ang, zeros(3, 1))
        # test set
        self.indiv.ang = Matrix([4, 5, 6])
        self.assertEqual(self.indiv.ang, Matrix([4, 5, 6]))
        with self.assertRaises(ShapeError):
            self.empty.ang = Matrix([3, 3])

    def test_equality(self):
        """Test __eq__() and __ne__()"""
        self.assertEqual(self.empty, Screw())
        self.assertNotEqual(self.indiv, Screw())
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
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(TestScrew)
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


