#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""Unit test for SymbolManager class."""


import unittest

from sympy import sympify, var, Matrix
from sympy.abc import A, B, C, X, Y, Z

from symoroutils import symbolmgr
from symoroutils import tools


class TestSymbolManager(unittest.TestCase):
    def setUp(self):
        self.symo = symbolmgr.SymbolManager()

    def test_get_max_coef(self):
        print("\n")
        expr1 = A*B*X + C**2 - X
        expr2 = Y*Z - B
        self.assertEqual(tools.get_max_coef(expr1*X + expr2, X), expr1)
        expr3 = -A**3*B**2*X**5*(X-Y)**7
        expr3x = -A**3*B**2*X**5*(-X-Y)**7
        expr3y = -A**3*B**2*X**5*(-X+Y)**7
        expr4 = B*X**2*(X-Y)**3
        self.assertEqual(tools.get_max_coef(expr3*expr4, expr4), expr3)
        self.assertEqual(tools.get_max_coef(expr3x, expr4), tools.ZERO)
        res = tools.get_max_coef(expr3y, expr4)*expr4-expr3y
        self.assertEqual(res.expand(), tools.ZERO)

    def test_name_extraction(self):
        print("\n")
        expr1 = sympify("C2*S3*R + S2*C3*R")
        self.assertEqual(tools.get_trig_couple_names(expr1), {'2', '3'})
        expr2 = sympify("CG2*S3*R + SG2*C1*R")
        self.assertEqual(tools.get_trig_couple_names(expr2), {'G2'})
        expr2 = sympify("CA2*SA3*R + SG2*C3*R")
        self.assertEqual(tools.get_trig_couple_names(expr2), set())
        expr3 = sympify("C2*S3*R + S1*C4*R")
        self.assertEqual(tools.get_trig_couple_names(expr3), set())

    def test_name_operations(self):
        print("\n")
        self.assertEqual(tools.reduce_str('12', '13'), ('2', '3'))
        self.assertEqual(tools.reduce_str('124', '123'), ('4', '3'))
        self.assertEqual(tools.reduce_str('124', '134'), ('2', '3'))
        self.assertEqual(tools.reduce_str('12', '124'), ('', '4'))
        self.assertEqual(tools.reduce_str('1G2', 'G24'), ('1', '4'))
        self.assertEqual(tools.reduce_str('1G2G4', '13G4'), ('G2', '3'))

    def test_try_opt(self):
        print("\n")
        e1 = A*(B-C)*X**2 + B*X**3 + A*(B-C)*Y**2 + B*X*Y**2
        e2 = X**2
        e3 = Y**2
        e4 = tools.ONE
        e5 = tools.ZERO
        self.assertEqual(self.symo.try_opt(e4, e5, e2, e3, e1), A*(B-C) + B*X)
        e6 = A*(B-C)*X**2 + B*X**3 - A*(B - C)*Y**2 - B*X*Y**2
        self.assertEqual(self.symo.try_opt(e4, e5, e2, e3, e6), e5)
        e7 = A*B
        self.assertEqual(self.symo.try_opt(e4, e7, e2, e3, e6),
                         e7*A*(B-C) + e7*B*X)
        self.assertEqual(self.symo.try_opt(e7, e4, e2, e3, e1),
                         e7*A*(B-C) + e7*B*X)

    def test_trig_simp(self):
        print("\n")
        e1 = sympify("S2**2 + C2**2")
        e1ans = sympify("1")
        self.assertEqual(self.symo.C2S2_simp(e1), e1ans)
        e1 = sympify("S1**2 + C2**2")
        self.assertEqual(self.symo.C2S2_simp(e1), e1)
        e1 = sympify("S2**3 + C2**2")
        self.assertEqual(self.symo.C2S2_simp(e1), e1)
        e1 = sympify("S2**2 + 2*C2**2")
        e1ans = sympify("C2**2 + 1")
        self.assertEqual(self.symo.C2S2_simp(e1), e1ans)
        e1 = sympify("S1**2 + S1**2*C1 + C1**2 + C1**3 + C1**4")
        e1ans = sympify("C1**4 + C1 + 1")
        self.assertEqual(self.symo.C2S2_simp(e1), e1ans)
        e2 = sympify("C1*S2 - C2*S1")
        e2ans = sympify("-S1m2")
        self.assertEqual(self.symo.CS12_simp(e2), e2ans)
        e2 = sympify("(C1*S2 - C2*S1)*(C1*S2 + C2*S1)")
        e2ans = sympify("-S1m2*S12")
        self.assertEqual(self.symo.CS12_simp(e2), e2ans)
        e2 = sympify("""C2*D3*S3m78 - C2m7*D8*S3 -
                     C3*D8*S2m7 - C3m78*D3*S2 + D2*S3""")
        e2ans = sympify("D2*S3 - D3*S278m3 - D8*S23m7")
        self.assertEqual(self.symo.CS12_simp(e2), e2ans)
        e2 = sympify("sin(g+th2)*sin(th3+th8)-cos(g+th2)*cos(th3+th8)")
        e2ans = sympify("-cos(g+th2+th3+th8)")
        self.assertEqual(self.symo.CS12_simp(e2), e2ans)
        e3 = sympify("""-a1*sin(th2+th1)*sin(th3)*cos(th1)-
                     a1*cos(th1)*cos(th2+th1)*cos(th3)""")
        e3ans = sympify("-a1*cos(th1)*cos(th1 + th2 - th3)")
        self.assertEqual(self.symo.CS12_simp(e3), e3ans)
        e4 = sympify("""C2*C3*C4**2*C5**2*C6**4*D3**2*RL4*S5 +
            2*C2*C3*C4**2*C5**2*C6**2*D3**2*RL4*S5*S6**2 +
            C2*C3*C4**2*C5**2*D3**2*RL4*S5*S6**4 +
            C2*C3*C4**2*C6**4*D3**2*RL4*S5**3 +
            2*C2*C3*C4**2*C6**2*D3**2*RL4*S5**3*S6**2 +
            C2*C3*C4**2*D3**2*RL4*S5**3*S6**4 +
            C2*C3*C5**2*C6**4*D3**2*RL4*S4**2*S5 +
            2*C2*C3*C5**2*C6**2*D3**2*RL4*S4**2*S5*S6**2 +
            C2*C3*C5**2*D3**2*RL4*S4**2*S5*S6**4 +
            C2*C3*C6**4*D3**2*RL4*S4**2*S5**3 +
            2*C2*C3*C6**2*D3**2*RL4*S4**2*S5**3*S6**2 +
            C2*C3*D3**2*RL4*S4**2*S5**3*S6**4 -
            C3*C4**2*C5**2*C6**4*D3*RL4**2*S23*S5 -
            2*C3*C4**2*C5**2*C6**2*D3*RL4**2*S23*S5*S6**2 -
            C3*C4**2*C5**2*D3*RL4**2*S23*S5*S6**4 -
            C3*C4**2*C6**4*D3*RL4**2*S23*S5**3 -
            2*C3*C4**2*C6**2*D3*RL4**2*S23*S5**3*S6**2 -
            C3*C4**2*D3*RL4**2*S23*S5**3*S6**4 -
            C3*C5**2*C6**4*D3*RL4**2*S23*S4**2*S5 -
            2*C3*C5**2*C6**2*D3*RL4**2*S23*S4**2*S5*S6**2 -
            C3*C5**2*D3*RL4**2*S23*S4**2*S5*S6**4 -
            C3*C6**4*D3*RL4**2*S23*S4**2*S5**3 -
            2*C3*C6**2*D3*RL4**2*S23*S4**2*S5**3*S6**2 -
            C3*D3*RL4**2*S23*S4**2*S5**3*S6**4""")
        e4ans = sympify("C3*D3*RL4*S5*(C2*D3 - RL4*S23)")
        self.assertEqual((self.symo.simp(e4)-e4ans).expand(), tools.ZERO)


def main():
    suite = unittest.TestLoader().loadTestsFromTestCase(
        TestSymbolManager
    )
    unittest.TextTestRunner(verbosity=2).run(suite)


if __name__ == '__main__':
    main()


