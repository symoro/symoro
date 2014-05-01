#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for TransformationMatrix class."""


import unittest

from pysymoro.transform import TransformationMatrix


class TestTransformationMatrix(unittest.TestCase):
    """Unit test for TransformationMatrix class."""
    pass


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


