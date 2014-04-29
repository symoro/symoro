#!/usr/bin/env python
# -*- coding: utf-8 -*-


import test_screw
import test_screw6
import test_geoparams
import test_dynparams


def main():
    """Main function."""
    test_screw.run_tests()
    test_screw6.run_tests()
    test_geoparams.run_tests()
    test_dynparams.run_tests()


if __name__ == '__main__':
    main()


