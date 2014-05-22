#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test module for GeoParams class."""


import unittest


from sympy import var, Matrix
from numpy import random, amax, matrix, eye, zeros


from pysymoro import invgeom
from symoroutils import samplerobots
from symoroutils import symbolmgr
from pysymoro import geometry


class TestIGM(unittest.TestCase):
    """Unit test for GeoParams class."""
    def setUp(self):
        self.symo = symbolmgr.SymbolManager()

    def test_igm_2r(self):
        print "######## test_igm ##########"
        robo = samplerobots.planar2r()
        nTm = Matrix(4, 4, 12 * [invgeom.EMPTY] + [0, 0, 0, 1])
        nTm[0, 3], nTm[1, 3] = var('p1, p2')
        invgeom._paul_solve(robo, self.symo, nTm, 0, robo.nf)
        self.symo.gen_func_string('IGM_gen', robo.q_vec,
                                  var('p1, p2'), syntax='matlab')
        igm_f = self.symo.gen_func('IGM_gen', robo.q_vec,
                                   var('p1, p2'))
        T = geometry.dgm(robo, self.symo, 0, robo.nf,
                         fast_form=True, trig_subs=True)
        f06 = self.symo.gen_func('DGM_generated1', (T[0, 3], T[1, 3]),
                                 robo.q_vec)
        for x in xrange(100):
            arg = random.normal(size=robo.nj)
            Ttest = f06(arg)
            solution = igm_f(Ttest)
            for q in solution:
                self.assertLess(amax(matrix(f06(q))-Ttest), 1e-12)

    def test_igm_rx90(self):
        print "######## test_igm ##########"
        robo = samplerobots.rx90()
        #robo.r[6] = var('R6')
        #robo.gamma[6] = var('G6')  # invgeom.T_GENERAL
        nTm = invgeom.T_GENERAL
        invgeom._paul_solve(robo, self.symo, nTm, 0, robo.nf)
        self.symo.gen_func_string('IGM_gen', robo.q_vec,
                                  invgeom.T_GENERAL, syntax='matlab')
        igm_f = self.symo.gen_func('IGM_gen', robo.q_vec,
                                   invgeom.T_GENERAL)
        T = geometry.dgm(robo, self.symo, 0, robo.nf,
                         fast_form=True, trig_subs=True)
        f06 = self.symo.gen_func('DGM_generated1', T, robo.q_vec)
        for x in xrange(100):
            arg = random.normal(size=robo.nj)
            Ttest = f06(arg)
            solution = igm_f(Ttest)
            for q in solution:
                self.assertLess(amax(matrix(f06(q))-Ttest), 1e-12)

    def test_loop(self):
        print "######## test_loop ##########"
        self.robo = samplerobots.sr400()
        invgeom.loop_solve(self.robo, self.symo)
        self.symo.gen_func_string('IGM_gen', self.robo.q_vec,
                                  self.robo.q_active, syntax='matlab')
        l_solver = self.symo.gen_func('IGM_gen', self.robo.q_vec,
                                      self.robo.q_active)
        T = geometry.dgm(self.robo, self.symo, 9, 10,
                         fast_form=True, trig_subs=True)
        t_loop = self.symo.gen_func('DGM_generated1', T, self.robo.q_vec)
        for x in xrange(10):
            arg = random.normal(size=6)
            solution = l_solver(arg)
            for q in solution:
                self.assertLess(amax(matrix(t_loop(q))-eye(4)), 1e-12)


def run_tests():
    """Load and run the unittests"""
    unit_suite = unittest.TestLoader().loadTestsFromTestCase(TestIGM)
    unittest.TextTestRunner(verbosity=2).run(unit_suite)


def main():
    """Main function."""
    run_tests()


if __name__ == '__main__':
    main()


