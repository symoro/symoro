#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for GeoParams class."""


import unittest

from sympy import var

from pysymoro.geoparams import GeoParams


class TestGeoParams(unittest.TestCase):
    """Unit test for GeoParams class."""
    def setUp(self):
        # setup an instance of GeoParams
        self.data_revol = GeoParams(22, {'sigma': 0})
        self.data_prism = GeoParams(20, {'sigma': 1})
        self.data_fixed = GeoParams(10, {'sigma': 2})
        # setup data to compare against
        self.frame = 33
        # setup values to compare against
        self.revol_theta = var('pi')
        self.prism_r = 5
        self.data_revol.theta = self.revol_theta
        self.data_prism.r = self.prism_r
        # setup params to compare against
        self.param_sigma = 0
        self.param_gamma = var('XX20')
        self.param_alpha = var('YY20')
        self.param_theta = var('ZZ20')
        self.params = {
            'sigma': self.param_sigma,
            'gamma': self.param_gamma,
            'alpha': self.param_alpha,
            'theta': self.param_theta
        }
        self.wrong_param = {'rand': 'some-value'}

    def test_init(self):
        """Test constructor"""
        # test instance type
        self.assertIsInstance(GeoParams(self.frame), GeoParams)

    def test_update_params(self):
        """Test update_params()"""
        # NOTE: the state of TransformationMatrix instance is not tested
        # here.
        # test raise AttributeError
        self.assertRaises(
            AttributeError, self.data_fixed.update_params,
            self.wrong_param
        )
        # test update values
        self.data_prism.update_params(self.params)
        self.assertEqual(self.data_prism.sigma, self.param_sigma)
        self.assertEqual(self.data_prism.gamma, self.param_gamma)
        self.assertEqual(self.data_prism.alpha, self.param_alpha)
        self.assertEqual(self.data_prism.theta, self.param_theta)

    def test_q(self):
        """Test get of q()"""
        self.assertEqual(self.data_revol.q, self.revol_theta)
        self.assertEqual(self.data_prism.q, self.prism_r)
        self.assertEqual(self.data_fixed.q, 0)


def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(
        TestGeoParams
    )
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


