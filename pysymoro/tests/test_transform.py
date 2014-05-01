#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for the functions in transform module."""


import unittest

from sympy import pi, var
from sympy import cos, sin
from sympy import Matrix

from pysymoro import transform


class TestTransform(unittest.TestCase):
    """Unit test for functions in transform module."""
    def setUp(self):
        # setup params for symbolic computation
        th1 = var('th1')
        l1 = var('L1')
        self.sym_gamma = 0
        self.sym_b = 0
        self.sym_alpha = 0
        self.sym_d = l1
        self.sym_theta = th1
        self.sym_r = 0
        self.sym_tmat = Matrix([
            [cos(th1), -sin(th1), 0, l1],
            [sin(th1), cos(th1), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        # setup params for numeric computation
        th1 = pi/2
        l1 = 1
        self.num_gamma = 0
        self.num_b = 0
        self.num_alpha = 0
        self.num_d = l1
        self.num_theta = th1
        self.num_r = 0
        self.num_tmat = Matrix([
            [cos(th1), -sin(th1), 0, l1],
            [sin(th1), cos(th1), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])

    def test_get_transformation_matrix_sym(self):
        """Test get_transformation_matrix() symbolically"""
        self.assertEqual(
            transform.get_transformation_matrix(
                self.sym_gamma, self.sym_b,
                self.sym_alpha, self.sym_d,
                self.sym_theta, self.sym_r
            ), self.sym_tmat
        )

    def test_get_transformation_matrix_num(self):
        """Test get_transformation_matrix() numerically"""
        self.assertEqual(
            transform.get_transformation_matrix(
                self.num_gamma, self.num_b,
                self.num_alpha, self.num_d,
                self.num_theta, self.num_r
            ), self.num_tmat
        )

    def test_get_transformation_matrix_sub(self):
        """Test get_transformation_matrix() by substitution"""
        t_in_sym = transform.get_transformation_matrix(
            self.sym_gamma, self.sym_b,
            self.sym_alpha, self.sym_d,
            self.sym_theta, self.sym_r
        )
        t_in_num = t_in_sym.subs({
            self.sym_d: self.num_d,
            self.sym_theta: self.num_theta
        })
        self.assertEqual(t_in_num, self.num_tmat)


def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(
        TestTransform
    )
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


