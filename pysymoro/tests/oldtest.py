#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
Unit tests for SYMORO modules
"""


import os
import unittest

from sympy import sympify, var, Matrix
from sympy.abc import A, B, C, X, Y, Z
from numpy import random, amax, matrix, eye, zeros

from pysymoro import robot
from pysymoro import geometry
from pysymoro.geometry import Transform as trns
from pysymoro import kinematics
from pysymoro import invgeom
from pysymoro import dynamics
from symoroutils import filemgr
from symoroutils import parfile
from symoroutils import samplerobots
from symoroutils import symbolmgr
from symoroutils import tools


class testMisc(unittest.TestCase):
    def test_robo_misc(self):
        print "######## test_robo_misc ##########"
        self.robo = samplerobots.sr400()
        q = list(var('th1:10'))
        self.assertEqual(self.robo.q_vec, q)
        self.assertEqual(self.robo.chain(6), [6, 5, 4, 3, 2, 1])
        self.assertEqual(self.robo.chain(6, 3), [6, 5, 4])
        self.assertEqual(self.robo.loop_chain(8, 9), [8, 9])
        self.assertEqual(self.robo.loop_chain(0, 6), [0, 1, 2, 3, 4, 5, 6])
        self.assertEqual(self.robo.loop_chain(6, 0), [6, 5, 4, 3, 2, 1, 0])
        self.assertEqual(self.robo.loop_chain(9, 10), [9, 8, 7, 1, 2, 3, 10])
        self.assertEqual(self.robo.loop_terminals, [(9, 10)])
        l1 = self.robo.get_geom_head()
        l2 = self.robo.get_dynam_head()
        l3 = self.robo.get_ext_dynam_head()
        for name in l1[1:] + l2[1:] + l3[1:]:
            for i in xrange(self.robo.NL):
                if name in tools.INT_KEYS:
                    self.assertEqual(self.robo.put_val(i, name, i), tools.OK)
                else:
                    v = var(name + str(i))
                    self.assertEqual(self.robo.put_val(i, name, v), tools.OK)
        for name in l3[1:]+l2[1:]+l1[1:]:
            for i in xrange(self.robo.NL):
                if name in tools.INT_KEYS:
                    self.assertEqual(self.robo.get_val(i, name), i)
                else:
                    v = var(name + str(i))
                    self.assertEqual(self.robo.get_val(i, name), v)


class testGeometry(unittest.TestCase):

    def setUp(self):
        self.symo = symbolmgr.SymbolManager()
        self.robo = samplerobots.rx90()

#    def test_misc(self):
#        self.assertEqual(self.robo.structure, tools.SIMPLE)
#        self.robo.ant[3] = 0
#        self.assertEqual(self.robo.type_of_structure, tools.TREE)
#        self.robo.ant[3] = 2
#        self.assertEqual(self.robo.type_of_structure, tools.SIMPLE)
#        robo2 = samplerobots.sr400()
#        self.assertEqual(robo2.type_of_structure, tools.CLOSED_LOOP)

    def test_dgm_rx90(self):
        print "######## test_dgm_rx90 ##########"
        T = geometry.dgm(self.robo, self.symo, 0, 6,
                         fast_form=True, trig_subs=True)
        f06 = self.symo.gen_func('DGM_generated1', T, self.robo.q_vec)
        T = geometry.dgm(self.robo, self.symo, 6, 0,
                         fast_form=True, trig_subs=True)
        f60 = self.symo.gen_func('DGM_generated2', T, self.robo.q_vec)
        for x in xrange(10):
            arg = random.normal(size=6)
            M = matrix(f06(arg))*matrix(f60(arg))-eye(4)
            self.assertLess(amax(M), 1e-12)
        t06 = matrix([[1, 0, 0, 1], [0, 1, 0, 0],
                      [0, 0, 1, 1], [0, 0, 0, 1]])
        self.assertLess(amax(matrix(f06(zeros(6)))-t06), 1e-12)
        T46 = geometry.dgm(self.robo, self.symo, 4, 6,
                           fast_form=False, trig_subs=True)
        C4, S4, C5, C6, S5, S6, RL4 = var("C4,S4,C5,C6,S5,S6,RL4")
        T_true46 = Matrix([[C5*C6, -C5*S6, -S5, 0], [S6, C6, 0, 0],
                           [S5*C6, -S5*S6, C5, 0], [0, 0, 0, 1]])
        self.assertEqual(T46, T_true46)
        T36 = geometry.dgm(self.robo, self.symo, 3, 6,
                           fast_form=False, trig_subs=True)
        T_true36 = Matrix([[C4*C5*C6-S4*S6, -C4*C5*S6-S4*C6, -C4*S5, 0],
                           [S5*C6, -S5*S6, C5, RL4],
                           [-S4*C5*C6-C4*S6, S4*C5*S6-C4*C6, S4*S5, 0],
                           [0, 0, 0, 1]])
        self.assertEqual(T36, T_true36)

    def test_dgm_sr400(self):
        print "######## test_dgm_sr400 ##########"
        self.robo = samplerobots.sr400()
        T = geometry.dgm(self.robo, self.symo, 0, 6,
                         fast_form=True, trig_subs=True)
        f06 = self.symo.gen_func('DGM_generated1', T, self.robo.q_vec)
        T = geometry.dgm(self.robo, self.symo, 6, 0,
                         fast_form=True, trig_subs=True)
        f60 = self.symo.gen_func('DGM_generated2', T, self.robo.q_vec)
        for x in xrange(10):
            arg = random.normal(size=9)
            M = matrix(f06(arg))*matrix(f60(arg))-eye(4)
            self.assertLess(amax(M), 1e-12)
        t06 = matrix([[1, 0, 0, 3], [0, -1, 0, 0],
                      [0, 0, -1, -1], [0, 0, 0, 1]])
        self.assertLess(amax(matrix(f06(zeros(9))) - t06), 1e-12)

    def test_igm(self):
        print "######## test_igm ##########"
        invgeom._paul_solve(self.robo, self.symo, invgeom.T_GENERAL, 0, 6)
        igm_f = self.symo.gen_func('IGM_gen', self.robo.q_vec,
                                   invgeom.T_GENERAL)
        T = geometry.dgm(self.robo, self.symo, 0, 6,
                         fast_form=True, trig_subs=True)
        f06 = self.symo.gen_func('DGM_generated1', T, self.robo.q_vec)
        for x in xrange(100):
            arg = random.normal(size=6)
            Ttest = f06(arg)
            solution = igm_f(Ttest)
            for q in solution:
                self.assertLess(amax(matrix(f06(q))-Ttest), 1e-12)

    def test_loop(self):
        print "######## test_loop ##########"
        self.robo = samplerobots.sr400()
        invgeom.loop_solve(self.robo, self.symo)
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


class testKinematics(unittest.TestCase):
    def setUp(self):
        self.symo = symbolmgr.SymbolManager()
        self.robo = samplerobots.rx90()

    def test_speeds(self):
        print 'Speeds and accelerations'
        kinematics.speeds_accelerations(self.robo)

        print 'Kinematic constraint equations'
        kinematics.kinematic_constraints(samplerobots.sr400())

    def test_jac(self):
        print "######## test_jac ##########"
        kinematics.jacobian(self.robo, 6, 3, 6)
        for j in xrange(1, 7):
            print "######## Jac validation through DGM ##########"
            #compute Jac
            J, l = kinematics._jac(self.robo, self.symo, j, 0, j)
            jacj = self.symo.gen_func('JacRX90', J, self.robo.q_vec)
            #compute DGM
            T = geometry.dgm(self.robo, self.symo, 0, j,
                             fast_form=True, trig_subs=True)
            T0j = self.symo.gen_func('DGM_generated1', T, self.robo.q_vec)
            for i in xrange(10):
                dq = random.normal(size=6, scale=1e-7)
                q = random.normal(size=6)
                dX = matrix(jacj(q)) * matrix(dq[:j]).T
                T = (matrix(T0j(q+dq)) - T0j(q))
                self.assertLess(amax(dX[:3] - trns.P(T)), 1e-12)

    def test_jac2(self):
        print "######## test_jac2 ##########"
        J, L = kinematics._jac(self.robo, self.symo, 6, 3, 3)
        jac63 = self.symo.gen_func('Jac1RX90', J, self.robo.q_vec)
        L63 = self.symo.gen_func('LRX90', L, self.robo.q_vec)
        J, L = kinematics._jac(self.robo, self.symo, 6, 3, 6)
        jac66 = self.symo.gen_func('Jac2RX90', J, self.robo.q_vec)
        for i in xrange(10):
            q = random.normal(size=6)
            j63 = matrix(jac63(q))
            l63 = matrix(L63(q))
            j66 = matrix(jac66(q))
            X = eye(6)
            X[:3, 3:] = l63
            self.assertLess(amax(j66 - X*j63), 1e-12)


class testDynamics(unittest.TestCase):
    def test_dynamics(self):
        robo = samplerobots.rx90()

        print 'Inverse dynamic model using Newton - Euler Algorith'
        dynamics.inverse_dynamic_NE(robo)

        print 'Direct dynamic model using Newton - Euler Algorith'
        dynamics.direct_dynamic_NE(robo)

        print 'Dynamic identification model using Newton - Euler Algorith'
        dynamics.dynamic_identification_NE(robo)

        print 'Inertia Matrix using composite links'
        dynamics.inertia_matrix(robo)

        print 'Coriolis torques using Newton - Euler Algorith'
        dynamics.pseudo_force_NE(robo)

        print 'Base parameters computation'
        dynamics.base_paremeters(robo)


if __name__ == '__main__':
    suite = unittest.TestSuite()
    suite.addTest(testMisc('test_robo_misc'))
    suite.addTest(testGeometry('test_dgm_rx90'))
    suite.addTest(testGeometry('test_dgm_sr400'))
    suite.addTest(testGeometry('test_igm'))
    suite.addTest(testGeometry('test_loop'))
    suite.addTest(testKinematics('test_jac'))
    suite.addTest(testKinematics('test_jac2'))
    unittest.TextTestRunner(verbosity=2).run(suite)


