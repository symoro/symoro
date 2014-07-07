#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module discovers all the unittest code in the folder and runs them.
"""


import inspect
import os
import unittest


def main():
    """Main function."""
    file_search_path = os.path.dirname(
        os.path.abspath(inspect.getfile(inspect.currentframe()))
    )
    all_tests = unittest.TestLoader().discover(
        file_search_path, pattern='test_*.py'
    )
    unittest.TextTestRunner(verbosity=2).run(all_tests)


if __name__ == '__main__':
    main()


