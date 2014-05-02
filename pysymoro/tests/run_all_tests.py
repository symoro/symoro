#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module discovers all the unittest code in the folder and runs them.
"""


import unittest


def main():
    """Main function."""
    all_tests = unittest.TestLoader().discover('tests', pattern='test_*.py')
    unittest.TextTestRunner(verbosity=2).run(all_tests)


if __name__ == '__main__':
    main()


