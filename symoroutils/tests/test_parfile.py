#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for functions in parfile.py file."""


import os
import unittest

from symoroutils import filemgr
from symoroutils import parfile
from symoroutils import samplerobots
from symoroutils import tools
from pysymoro import robot


class TestParfile(unittest.TestCase):
    def setUp(self):
        self.orig_robo = samplerobots.rx90()

    def test_readwrite(self):
        parfile.writepar(self.orig_robo)
        fname = filemgr.get_clean_name(self.orig_robo.name) + ".par"
        file_path = os.path.join(self.orig_robo.directory, fname)
        new_robo, flag = parfile.readpar(self.orig_robo.name, file_path)
        self.assertEqual(flag, tools.OK)
        l1 = self.orig_robo.get_geom_head()
        l2 = self.orig_robo.get_dynam_head()
        l3 = self.orig_robo.get_ext_dynam_head()
        for name in l3[1:]+l2[1:]+l1[1:]:
            for i in xrange(1, self.orig_robo.NL):
                self.assertEqual(self.orig_robo.get_val(i, name),
                                 new_robo.get_val(i, name))


def main():
    # load and run the doctests
#    doc_suite = doctest.DocTestSuite(parfile)
#    unittest.TextTestRunner(verbosity=2).run(doc_suite)
    # load and run the unittests
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(TestParfile)
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


if __name__ == '__main__':
    main()


