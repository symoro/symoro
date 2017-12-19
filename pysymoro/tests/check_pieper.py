#!/usr/bin/env python
# -*- coding: utf-8 -*-

import unittest

from sympy import var, Matrix
from numpy import random, amax, matrix, eye, zeros

from pysymoro import geometry, pieper, invdata
from symoroutils import samplerobots, symbolmgr

class TestIGM(unittest.TestCase):
    """Unit test for GeoParams class."""
    def setUp(self):
        self.symo = symbolmgr.SymbolManager()

    def test_igm_rx90(self):
        robo = samplerobots.rx90()
        nTm = invdata.T_GENERAL
        self.symo.write_line("\r\n*******************PIEPER METHOD STARTS********************** \r\n")
        pieper._pieper_solve(robo, self.symo, nTm, 0, robo.nf)
##        self.symo.gen_func_string('IGM_gen', robo.q_vec,
##                                  invgeom.T_GENERAL, syntax='matlab')
        igm_f = self.symo.gen_func('IGM_gen', robo.q_vec,
                                   invdata.T_GENERAL)
        self.symo.write_line("\r\n*******************PIEPER METHOD FINISHED********************")
##        T = geometry.dgm(robo, self.symo, 0, robo.nf,
##                         fast_form=True, trig_subs=True)
##        f06 = self.symo.gen_func('DGM_generated1', T, robo.q_vec)
##        for x in range(100):
##            arg = random.normal(size=robo.nj)
##            Ttest = f06(arg)
##            solution = igm_f(Ttest)

def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(TestIGM)
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()
