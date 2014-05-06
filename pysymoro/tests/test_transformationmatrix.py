#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for TransformationMatrix class."""


import unittest

from sympy import eye, var, zeros
from sympy import cos, sin
from sympy import Matrix

from pysymoro.screw6 import Screw6
from symoroutils import tools

from pysymoro.transform import TransformationMatrix


class TestTransformationMatrix(unittest.TestCase):
    """Unit test for TransformationMatrix class."""
    def setUp(self):
        # parameter values
        th1 = var('th1')
        l1 = var('L1')
        self.gamma = 0
        self.b = 0
        self.alpha = 0
        self.d = l1
        self.theta = th1
        self.r = 0
        self.t_val = Matrix([
            [cos(th1), -sin(th1), 0, l1],
            [sin(th1), cos(th1), 0, 0],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        self.rot_val = Matrix([
            [cos(th1), -sin(th1), 0],
            [sin(th1), cos(th1), 0],
            [0, 0, 1]
        ])
        self.trans_val = Matrix([l1, 0, 0])
        self.t_inv = Matrix([
            [cos(th1), sin(th1), 0, -l1*cos(th1)],
            [-sin(th1), cos(th1), 0, l1*sin(th1)],
            [0, 0, 1, 0],
            [0, 0, 0, 1]
        ])
        self.inv_rot = Matrix([
            [cos(th1), sin(th1), 0],
            [-sin(th1), cos(th1), 0],
            [0, 0, 1]
        ])
        self.inv_trans = Matrix([-l1*cos(th1), l1*sin(th1), 0])
        self.sji = Screw6(
            tl=self.rot_val, tr=tools.skew(self.trans_val),
            bl=zeros(3, 3), br=self.rot_val
        ).val
        self.sij = Screw6(
            tl=self.inv_rot,
            tr=-(self.inv_rot * tools.skew(self.trans_val)),
            bl=zeros(3, 3), br=self.inv_rot
        ).val
        # params dict
        self.params = {
            'gamma': self.gamma,
            'b': self.b,
            'alpha': self.alpha,
            'd': self.d,
            'theta': self.theta,
            'r': self.r
        }
        self.sigmu_param = {'sigma': 2, 'mu': 0}
        self.wrong_param = {'rand': 'some-value'}
        # frames
        self.frame_i = 5
        self.frame_j = 6
        # setup instances for use
        self.t_empty = TransformationMatrix(i=0, j=0)
        self.t_data = TransformationMatrix(
            i=self.frame_i, j=self.frame_j, params=self.params
        )

    def test_init(self):
        """Test contructor"""
        params = self.params
        # test raise NotImplementedError
        self.assertRaises(NotImplementedError, TransformationMatrix)
        # test raise AttributeError
        self.assertRaises(
            AttributeError, TransformationMatrix, params=params
        )
        # test different forms of the constructor
        self.assertIsInstance(
            TransformationMatrix(i=0, j=0), TransformationMatrix
        )
        self.assertIsInstance(
            TransformationMatrix(i=0, j=0, params=params),
            TransformationMatrix
        )
        params['frame'] = 1
        params['ant'] = 0
        self.assertIsInstance(
            TransformationMatrix(params=params), TransformationMatrix
        )

    def test_val(self):
        """Test get of val()"""
        self.assertEqual(self.t_empty.val, eye(4))
        self.assertEqual(self.t_data.val, self.t_val)

    def test_rot(self):
        """Test get of rot()"""
        self.assertEqual(self.t_empty.rot, eye(3))
        self.assertEqual(self.t_data.rot, self.rot_val)

    def test_trans(self):
        """Test get of trans()"""
        self.assertEqual(self.t_empty.trans, zeros(3, 1))
        self.assertEqual(self.t_data.trans, self.trans_val)

    def test_inv(self):
        """Test get of inv()"""
        self.assertEqual(self.t_empty.inv, eye(4))
        self.assertEqual(self.t_data.inv, self.t_inv)

    def test_inv_rot(self):
        """Test get of inv_rot()"""
        self.assertEqual(self.t_empty.inv_rot, eye(3))
        self.assertEqual(self.t_data.inv_rot, self.inv_rot)

    def test_inv_trans(self):
        """Test get of inv_trans()"""
        self.assertEqual(self.t_empty.inv_trans, zeros(3, 1))
        self.assertEqual(self.t_data.inv_trans, self.inv_trans)

    def test_s_j_wrt_i(self):
        """Test get of s_j_wrt_i()"""
        self.assertEqual(self.t_empty.s_j_wrt_i, eye(6))
        self.assertEqual(self.t_data.s_j_wrt_i, self.sji)

    def test_s_i_wrt_j(self):
        """Test get of s_i_wrt_j()"""
        self.assertEqual(self.t_empty.s_i_wrt_j, eye(6))
        self.assertEqual(self.t_data.s_i_wrt_j, self.sij)

    def test_update(self):
        """Test update()"""
        # test raise AttributeError
        self.assertRaises(
            AttributeError, self.t_empty.update, self.wrong_param
        )
        # test sigma, mu
        self.t_empty.update(self.sigmu_param)
        self.assertEqual(self.t_empty.val, eye(4))
        # test with correct params
        self.t_empty.update(self.params)
        self.assertEqual(self.t_empty.val, self.t_val)


def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(
        TestTransformationMatrix
    )
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


