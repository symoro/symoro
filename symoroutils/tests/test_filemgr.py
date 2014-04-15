#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for functions in filemgr.py file."""


import doctest
import os
import shutil
import unittest
from collections import namedtuple

from symoroutils import filemgr


class TestFileMgr(unittest.TestCase):
    """Unit test for functions in filemgr.py file"""
    def setUp(self):
        # strings that will be used by different test cases
        self.home_folder = os.path.expanduser("~")
        self.robots_folder = "symoro-robots"
        self.base_path = os.path.join(
            self.home_folder, self.robots_folder
        )
        self.rob1 = "File Test Robot"
        self.rob2 = "FileTestRobot"
        # tmp robot created using namedtuple just for the purpose of
        # test case
        RobotTmp = namedtuple('RobotTmp', ['name', 'directory'])
        self.tmp_robot = RobotTmp('Tmp Robot', self.base_path)

    def test_get_base_path(self):
        self.assertEqual(filemgr.get_base_path(), self.base_path)

    def test_get_folder_path(self):
        rob1_clean = filemgr.get_clean_name(self.rob1)
        rob2_clean = filemgr.get_clean_name(self.rob2)
        # scenario 1
        self.assertEqual(
            filemgr.get_folder_path(self.rob1),
            os.path.join(self.base_path, rob1_clean)
        )
        # scenario 2
        self.assertEqual(
            filemgr.get_folder_path(self.rob2),
            os.path.join(self.base_path, rob2_clean)
        )

    def test_make_file_path(self):
        robot_name_clean = filemgr.get_clean_name(self.tmp_robot.name)
        # scenario 1
        par_checker = os.path.join(
            self.tmp_robot.directory,
            '%s.par' % robot_name_clean
        )
        par_format = filemgr.make_file_path(self.tmp_robot)
        self.assertEqual(par_format, par_checker)
        # scenario 2
        trm_checker = os.path.join(
            self.tmp_robot.directory,
            '%s_%s.txt' % (robot_name_clean, 'trm')
        )
        trm_format = filemgr.make_file_path(self.tmp_robot, 'trm')
        self.assertEqual(trm_format, trm_checker)

    def tearDown(self):
        # delete the temp folders and files created
        rob1_path = filemgr.get_folder_path(self.rob1)
        rob2_path = filemgr.get_folder_path(self.rob2)
        if os.path.exists(rob1_path):
            shutil.rmtree(rob1_path)
        if os.path.exists(rob2_path):
            shutil.rmtree(rob2_path)


def main():
    # load and run the doctests
    doc_suite = doctest.DocTestSuite(filemgr)
    unittest.TextTestRunner(verbosity=2).run(doc_suite)
    # load and run the unittests
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(TestFileMgr)
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


if __name__ == '__main__':
    main()


