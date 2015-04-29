# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package provides symbolic
solutions for inverse geompetric problem fro robots with a simple geometry.
"""

from heapq import heapify, heappop
from itertools import combinations

from sympy import var, sin, cos, atan2, sqrt, pi
from sympy import Matrix, Symbol, Expr
from pysymoro.geometry import TransformChain
from symoroutils.tools import ZERO, ONE
from symoroutils import symbolmgr

EMPTY = var("EMPTY")
T_GENERAL = Matrix([var("s1,n1,a1,p1"), var("s2,n2,a2,p2"),
                    var("s3,n3,a3,p3"), [0, 0, 0, 1]])

#TODO: revamp the notation, remove trigsimp
theta1, theta2, rho = var("theta1, theta2, rho")

C1 = cos(theta1)
C2 = cos(theta2)
S1 = sin(theta1)
S2 = sin(theta2)
C12 = cos(theta1 + theta2)
S12 = sin(theta1 + theta2)


def trigsimp(A):
    return A


def replace_EMPTY(T, tr_list):
    T_sym = tr_list.to_matrix(simplify=True)
    for e1 in xrange(4):
        for e2 in xrange(4):
            if T[e1, e2].has(EMPTY):
                T[e1, e2] = T_sym[e1, e2]


def loop_solve(robo, symo, know=None):
    # TODO: rewrite; Add parallelogram detection
    q_vec = [robo.get_q(i) for i in xrange(robo.NF)]
    loops = []
    if know is None:
        know = robo.q_active
        # set(q for i, q in enumerate(q_vec) if robo.mu[i] == 1)
    for i, j in robo.loop_terminals:
        chain = robo.loop_chain(i, j)
        know_ij = set(q_vec[i] for i in chain if q_vec[i] in know)
        unknow_ij = set(q_vec[i] for i in chain if q_vec[i] not in know)
        loops.append([len(unknow_ij), i, j, know_ij, unknow_ij])
    while loops:
        heapify(loops)
        loop = heappop(loops)
        res_know = paul_solve(robo, symo, eye(4), *loop[1:4])
        for l in loops:
            found = l[4] & res_know
            l[3] |= found
            l[4] -= found
            l[0] = len(l[4])


def igm_paul(robo, T_ref, n):
    if isinstance(T_ref, list):
        T_ref = Matrix(4, 4, T_ref)
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'igm')
    symo.write_params_table(robo, 'Inverse Geometric Model for frame %s' % n)
    paul_solve(robo, symo, T_ref, 0, n)
    symo.file_close()
    return symo


def paul_solve(robo, symo, nTm, n, m, known_vars=None):
    chain = robo.loop_chain(m, n)
    th_all = set()
    r_all = set()
    # Create the set of all knowns symbols
    for i in chain:
        if i >= 0:
            if robo.sigma[i] == 0 and isinstance(robo.theta[i], Expr):
                th_all.add(robo.theta[i])
            if robo.sigma[i] == 1 and isinstance(robo.r[i], Expr):
                r_all.add(robo.r[i])
    solver = Solver(symo, th_all, r_all)
    prev_unknowns_size = solver.unknowns_size
    while not solver.is_solved:
        tr_list = TransformChain.create(robo, n, m)
        Td = nTm.copy()
        replace_EMPTY(Td, tr_list)
        #bring all the terminal constant transforms on the right side
        tr_const = tr_list.separate_const_right(solver.unknowns)
        tr_post = tr_const.inverse()
        Td = Td * tr_post.to_matrix(simplify=True)
        tr_pre = TransformChain()
        while tr_list.data:
            print "ITERATION"
            tr_const = tr_list.separate_const_left(solver.unknowns)
            assert isinstance(tr_pre, TransformChain)
            tr_pre = tr_const.inverse() + tr_pre
            assert isinstance(tr_pre, TransformChain)
            Tleft = tr_list.to_matrix(simplify=True)
            Tright = tr_pre.to_matrix(simplify=True)*Td
            M_eq = Tleft - Tright
            tr = tr_list.pop(0)
            tr_pre.insert(0, tr.inv())
           # if tr.val.atoms(Symbol) - knowns:
            print "### EQUATIONS ###"
            for eq in list(M_eq[0:3, 0:4]):
                print eq
                print '---'
            solver.solve(list(M_eq[0:3, 0:4]))
        if prev_unknowns_size == solver.unknowns_size:
            break


class Pattern:
    """
    Contains a description of a certain solvable equation
    """
    def __init__(self, is_constant=None):
        #list of variable terms
        self.var_lists = []
        #constraints for double equation system
        self.equal_constr = []
        self.inverse_constr = []
        #for terms that are on the right hand side
        self.negative_coeffs = []
        self.is_constant = is_constant

    def match(self, *equation_list, **kwvar):
        if len(equation_list) != len(self.var_lists):
            return None
        coeff_list = []
        #TODO: is it possible to do the equation processing
        # just once. Right now it is here and in prepare_equations
        for i, equation in enumerate(equation_list):
            if "r" in kwvar:
                equation = equation.subs(kwvar["r"], rho)
            if "th1" in kwvar:
                equation = equation.subs(kwvar["th1"], theta1)
            if "th2" in kwvar:
                equation = equation.subs(kwvar["th2"], theta2)
            coeffs = self.exctract_coeffs(equation, i)
            if not coeffs:
                return None
            coeff_list.append(coeffs)
        if len(coeff_list) == 1:
            return coeff_list[0]
        elif self.check_constr(coeff_list[0], coeff_list[1]):
            skip_coeffs = set(self.equal_constr) | set(self.inverse_constr)
            for i, coeff in enumerate(coeff_list[1]):
                if i not in skip_coeffs:
                    coeff_list[0].append(coeff)
            return coeff_list[0]
        else:
            return None

    def exctract_coeffs(self, equation, var_idx):
        coeffs = []
        rest = equation
        for var in self.var_lists[var_idx]:
            rest, term = rest.as_independent(var, as_Add=True)
            if term == ZERO:
                coeffs.append(ZERO)
                continue
            c = term.as_coefficient(var)
            if c is None or not self.is_constant(c):
                return None
            coeffs.append(c)
        #constant term
        if not self.is_constant(rest):
            return None
        coeffs.append(-rest)
        for i in self.negative_coeffs:
            coeffs[i] = -coeffs[i]
        return coeffs

    def check_constr(self, coeffs1, coeffs2):
        positive = True
        for i in self.equal_constr:
            if coeffs1[i] != coeffs2[i]:
                positive = False
        for i in self.inverse_constr:
            if coeffs1[i] != -coeffs2[i]:
                positive = False
        if positive:
            return True
        for i in self.equal_constr:
            if coeffs1[i] != -coeffs2[i]:
                return False
        for i in self.inverse_constr:
            if coeffs1[i] != coeffs2[i]:
                return False
        for i in xrange(len(coeffs2)):
            coeffs2[i] = -coeffs2[i]
        return True


class Solver:
    """
    Solves several types of equations and systems of equations
    solve(...) performs search across given list of equations
    in order to find known type and solve it.
    """
    def __init__(self, symo, unknown_th, unknown_r=set()):
        self.symo = symo
        self.unknown_th = set(unknown_th)
        self.unknown_r = set(unknown_r)
        #initialize patterns
        self.patt = {}
        for i in (1, 2, 3, 4, 5, 6, 7, 8):
            self.patt[i] = Pattern(self.is_constant)
        self.patt[1].var_lists = [[rho]]
        self.patt[2].var_lists = [[S1, C1]]
        self.patt[3].var_lists = [[S1, C1], [S1, C1]]
        self.patt[4].var_lists = [[S1*rho], [C1*rho]]
        self.patt[5].var_lists = [[S1, rho], [C1, rho]]
        self.patt[5].negative_coeffs = [1]
        self.patt[6].var_lists = [[S1, C1, C2, S2], [C1, S1, S2, C2]]
        self.patt[6].negative_coeffs = [2, 3]
        self.patt[6].equal_constr = [0, 2]
        self.patt[6].inverse_constr = [1, 3]
        self.patt[7].var_lists = [[C1, S1, C2, S2], [S1, C1, S2, C2]]
        self.patt[7].negative_coeffs = [2, 3]
        self.patt[7].equal_constr = [0, 2]
        self.patt[7].inverse_constr = [1, 3]
        self.patt[8].var_lists = [[C1, C12], [S1, S12]]
        self.patt[8].equal_constr = [0, 1]

    @property
    def is_solved(self):
        if self.unknown_th or self.unknown_r:
            return False
        else:
            return True

    @property
    def unknowns_size(self):
        return len(self.unknown_th) + len(self.unknown_r)

    @property
    def unknowns(self):
        return self.unknown_th | self.unknown_r

    def prepare_equations(self, equation_list):
        """
        TO WRITE A DESCRIPTION
        """
        equation_candidates = []
        for equation in equation_list:
            ths, rs = self.retrieve_unknowns(equation)
            if len(ths | rs) > 2 or len(ths | rs) == 0 or len(rs) > 1:
                continue
            collect_terms = []
            if len(rs) == 1:
                q, = rs
                collect_terms.append(q)
            if len(ths) == 2:
                q1, q2 = ths
                collect_terms.append(cos(q1))
                collect_terms.append(cos(q2))
                collect_terms.append(sin(q1))
                collect_terms.append(sin(q2))
                collect_terms.append(cos(q1+q2))
                collect_terms.append(sin(q1+q2))
            elif len(ths) == 1:
                q, = ths
                collect_terms.append(cos(q))
                collect_terms.append(sin(q))
            #equation = equation.expand()
            equation = equation.collect(collect_terms)
            equation_candidates.append((ths, rs, equation))
        return equation_candidates

    def solve(self, equation_list):
        """
        Looks for euqations according to the known patterns
        and calls functions to write the solution.
        """
        equation_candidates = self.prepare_equations(equation_list)
        #look for double equations
        for eq_pack1, eq_pack2 in combinations(equation_candidates, 2):
            ths = eq_pack1[0] | eq_pack2[0]
            rs = eq_pack1[1] | eq_pack2[1]
            if len(ths | rs) > 2:
                continue
            if self.solve_double(eq_pack1[2], eq_pack2[2], ths, rs):
                return True
        #look for a single solvable equation
        for eq_pack in equation_candidates:
            ths = eq_pack[0]
            rs = eq_pack[1]
            if len(ths | rs) > 1:
                continue
            if self.solve_mono(eq_pack[2], ths, rs):
                return True
        return False

    def solve_double(self, eq1, eq2, ths, rs):
        """
        Tries to classify the equation system and call a fucntion that writes
        a solution.
        """
        if len(ths) == 1 and len(rs) == 0:
            th, = ths
            coeffs = self.patt[3].match(eq1, eq2, th1=th)
            if coeffs:
                _solve_type_3(self.symo, coeffs, th)
                self.unknown_th.remove(th)
                print "TYPE 3"
                return True
        elif len(ths) == 1 and len(rs) == 1:
            th, = ths
            r, = rs
            coeffs = self.patt[4].match(eq1, eq2, th1=th, r=r)
            if not coeffs:  # try the other way around
                coeffs = self.patt[4].match(eq2, eq1, th1=th, r=r)
            if coeffs:
                _solve_type_4(self.symo, coeffs, th, r)
                self.unknown_th.remove(th)
                self.unknown_r.remove(r)
                print "TYPE 4"
                return True
            coeffs = self.patt[5].match(eq1, eq2, th1=th, r=r)
            if not coeffs:  # try the other way around
                coeffs = self.patt[5].match(eq2, eq1, th1=th, r=r)
            if coeffs:
                _solve_type_5(self.symo, coeffs, th, r)
                self.unknown_th.remove(th)
                self.unknown_r.remove(r)
                print "TYPE 5"
                return True
        elif len(ths) == 2 and len(rs) == 0:
            th1, th2 = ths
            coeffs = self.patt[6].match(eq1, eq2, th1=th1, th2=th2)
            if coeffs:
                _solve_type_6(self.symo, coeffs, th1, th2)
                self.unknown_th.remove(th1)
                self.unknown_th.remove(th2)
                print "TYPE 6"
                return True
            coeffs = self.patt[7].match(eq1, eq2, th1=th1, th2=th2)
            if coeffs:
                _solve_type_7(self.symo, coeffs, th1, th2)
                self.unknown_th.remove(th1)
                self.unknown_th.remove(th2)
                print "TYPE 7"
                return True
            coeffs = self.patt[8].match(eq1, eq2, th1=th1, th2=th2)
            if not coeffs:  # try the other way around
                th1, th2 = th2, th1
                coeffs = self.patt[8].match(eq1, eq2, th1=th1, th2=th2)
            if coeffs:
                _solve_type_8(self.symo, coeffs, th1, th2)
                self.unknown_th.remove(th1)
                self.unknown_th.remove(th2)
                print "TYPE 8"
                return True

    def solve_mono(self, eq, ths, rs):
        """
        Tries to classify the equation and call a fucntion that writes
        a solution.
        """
        if len(ths) == 0 and len(rs) == 1:
            r, = rs
            coeffs = self.patt[1].match(eq, r=r)
            if coeffs:
                _solve_type_1(self.symo, coeffs, r)
                self.unknown_r.remove(r)
                print "TYPE 1"
                return True
        elif len(ths) == 1 and len(rs) == 0:
            th, = ths
            coeffs = self.patt[2].match(eq, th1=th)
            if coeffs:
                _solve_type_2(self.symo, coeffs, th)
                self.unknown_th.remove(th)
                print "TYPE 2"
                return True
        return False

    def retrieve_unknowns(self, equation):
        """
        Returns unknown thetas and rs.
        """
        thetas = equation.atoms(Symbol) & self.unknown_th
        rs = equation.atoms(Symbol) & self.unknown_r
        return thetas, rs

    def is_constant(self, equation):
        if not isinstance(equation, Expr):
            return True
        eq_syms = equation.atoms(Symbol)
        for var_sym in (theta1, theta2, rho):
            if var_sym in eq_syms:
                return False
        th, r = self.retrieve_unknowns(equation)
        if th or r:
            return False
        else:
            return True


def _solve_type_1(symo, coeffs, r):
    """Solution for the equation:
    X*r = Y
    """
    X, Y = coeffs
    symo.write_line("# X*r = Y".format(r))
    X = symo.replace(trigsimp(X), 'X', r)
    Y = symo.replace(trigsimp(Y), 'Y', r)
    symo.add_to_dict(r, Y/X)


def _solve_type_2(symo, coeffs, th):
    """Solution for the equation:
    X*S + Y*C = Z
    """
    X, Y, Z = coeffs
    symo.write_line("# X*sin({0}) + Y*cos({0}) = Z".format(th))
    X = symo.replace(trigsimp(X), 'X', th)
    Y = symo.replace(trigsimp(Y), 'Y', th)
    Z = symo.replace(trigsimp(Z), 'Z', th)
    YPS = var('YPS%s' % th)
    if X == ZERO and Y != ZERO:
        C = symo.replace(Z/Y, 'C', th)
        symo.add_to_dict(YPS, (ONE, - ONE))
        symo.add_to_dict(th, atan2(YPS*sqrt(1-C**2), C))
    elif X != ZERO and Y == ZERO:
        S = symo.replace(Z/X, 'S', th)
        symo.add_to_dict(YPS, (ONE, - ONE))
        symo.add_to_dict(th, atan2(S, YPS*sqrt(1-S**2)))
    elif Z == ZERO:
        symo.add_to_dict(YPS, (ONE, ZERO))
        symo.add_to_dict(th, atan2(-Y, X) + YPS*pi)
    else:
        B = symo.replace(X**2 + Y**2, 'B', th)
        D = symo.replace(B - Z**2, 'D', th)
        symo.add_to_dict(YPS, (ONE, - ONE))
        S = symo.replace((X*Z + YPS * Y * sqrt(D))/B, 'S', th)
        C = symo.replace((Y*Z - YPS * X * sqrt(D))/B, 'C', th)
        symo.add_to_dict(th, atan2(S, C))


def _solve_type_3(symo, coeffs, th):
    """Solution for the system:
    X1*S + Y1*C = Z1
    X2*S + Y2*C = Z2
    """
    X1, Y1, Z1, X2, Y2, Z2 = coeffs
    symo.write_line("# X1*sin({0}) + Y1*cos({0}) = Z1".format(th))
    symo.write_line("# X2*sin({0}) + Y2*cos({0}) = Z2".format(th))
    X1 = symo.replace(trigsimp(X1), 'X1', th)
    Y1 = symo.replace(trigsimp(Y1), 'Y1', th)
    Z1 = symo.replace(trigsimp(Z1), 'Z1', th)
    X2 = symo.replace(trigsimp(X2), 'X2', th)
    Y2 = symo.replace(trigsimp(Y2), 'Y2', th)
    Z2 = symo.replace(trigsimp(Z2), 'Z2', th)
    if X1 == ZERO and Y2 == ZERO:
        symo.add_to_dict(th, atan2(Z2/X2, Z1/Y1))
    elif X2 == ZERO and Y1 == ZERO:
        symo.add_to_dict(th, atan2(Z1/X1, Z2/Y2))
    else:
        D = symo.replace(X1*Y2-X2*Y1, 'D', th)
        C = symo.replace((Z2*X1 - Z1*X2)/D, 'C', th)
        S = symo.replace((Z1*Y2 - Z2*Y1)/D, 'S', th)
        symo.add_to_dict(th, atan2(S, C))


def _solve_type_4(symo, coeffs, th, r):
    """Solution for the system:
    X1*S*r = Y1
    X2*C*r = Y2
    """
    X1, Y1, X2, Y2 = coeffs
    symo.write_line("# X1*sin({0})*{1} = Y1".format(th, r))
    symo.write_line("# X2*cos({0})*{1} = Y2".format(th, r))
    X1 = symo.replace(trigsimp(X1), 'X1', th)
    Y1 = symo.replace(trigsimp(Y1), 'Y1', th)
    X2 = symo.replace(trigsimp(X2), 'X2', th)
    Y2 = symo.replace(trigsimp(Y2), 'Y2', th)
    YPS = var('YPS%s' % r)
    symo.add_to_dict(YPS, (ONE, - ONE))
    symo.add_to_dict(r, YPS*sqrt((Y1/X1)**2 + (Y2/X2)**2))
    symo.add_to_dict(th, atan2(Y1/(X1*r), Y2/(X2*r)))


def _solve_type_5(symo, coeffs, th, r):
    """Solution for the system:
    X1*S = Y1*r + Z1
    X2*C = Y2*r + Z2
    """
    X1, Y1, Z1, X2, Y2, Z2 = coeffs
    symo.write_line("# X1*sin({0}) = Y1*{1} + Z1".format(th, r))
    symo.write_line("# X2*cos({0}) = Y2*{1} + Z2".format(th, r))
    X1 = symo.replace(trigsimp(X1), 'X1', th)
    Y1 = symo.replace(trigsimp(Y1), 'Y1', th)
    Z1 = symo.replace(trigsimp(Z1), 'Z1', th)
    X2 = symo.replace(trigsimp(X2), 'X2', th)
    Y2 = symo.replace(trigsimp(Y2), 'Y2', th)
    Z2 = symo.replace(trigsimp(Z2), 'Z2', th)
    V1 = symo.replace(Z1/X1, 'V1', r)
    W1 = symo.replace(Y1/X1, 'W1', r)
    V2 = symo.replace(Z2/X2, 'V2', r)
    W2 = symo.replace(Y2/X2, 'W2', r)
    _solve_square(W1**2 + W2**2, 2*(V1*W1 + V2*W2), V1**2 + V2**2, r)
    _solve_type_3(symo, [X1, ZERO, Y1*r + Z1, ZERO, X2, Y2*r + Z2], th)


#TODO: compress with type_7
def _solve_type_6(symo, coeffs, th_1, th_2):
    """Solution for the system:
    W*S1 + V*C1 = X*C2 + Y*S2 + Z1
    W*C1 - V*S1 = X*S2 - Y*C2 + Z2
    """
    W, V, X, Y, Z1, Z2 = coeffs
    s = "# V*cos({0}) + W*sin({0}) = X*cos({1}) + Y*sin({1}) + Z1"
    symo.write_line(s.format(th_1, th_2))
    s = "# V*sin({0}) - W*cos({0}) = X*sin({1}) - Y*cos({1}) + Z2"
    symo.write_line(s.format(th_1, th_2))
    V = symo.replace(trigsimp(V), 'V', th_2)
    W = symo.replace(trigsimp(W), 'W', th_2)
    X = symo.replace(trigsimp(X), 'X', th_2)
    Y = symo.replace(trigsimp(Y), 'Y', th_2)
    Z1 = symo.replace(trigsimp(Z1), 'Z1', th_2)
    Z2 = symo.replace(trigsimp(Z2), 'Z2', th_2)
    B1 = symo.replace(2*(Z1*Y + Z2*X), 'B1', th_2)
    B2 = symo.replace(2*(Z1*X - Z2*Y), 'B2', th_2)
    B3 = symo.replace(V**2 + W**2 - X**2 - Y**2 - Z1**2 - Z2**2, 'B3', th_2)
    _solve_type_2(symo, [B1, B2, B3], th_2)
    Zi1 = symo.replace(X*cos(th_2) + Y*sin(th_2) + Z1, 'Zi1', th_1)
    Zi2 = symo.replace(X*sin(th_2) - Y*cos(th_2) + Z2, 'Zi2', th_1)
    _solve_type_3(symo, [W, V, Zi1, -V, W, Zi2], th_1)


def _solve_type_7(symo, coeffs, th_1, th_2):
    """Solution for the system:
    V*Cj + W*Sj = X*Ci + Y*Si + Z1
    V*Sj - W*Cj = X*Si - Y*Ci + Z2
    """
    V, W, X, Y, Z1, Z2 = coeffs
    s = "# V*cos({0}) + W*sin({0}) = X*cos({1}) + Y*sin({1}) + Z1"
    symo.write_line(s.format(th_1, th_2))
    s = "# V*sin({0}) - W*cos({0}) = X*sin({1}) - Y*cos({1}) + Z2"
    symo.write_line(s.format(th_1, th_2))
    V = symo.replace(trigsimp(V), 'V', th_2)
    W = symo.replace(trigsimp(W), 'W', th_2)
    X = symo.replace(trigsimp(X), 'X', th_2)
    Y = symo.replace(trigsimp(Y), 'Y', th_2)
    Z1 = symo.replace(trigsimp(Z1), 'Z1', th_2)
    Z2 = symo.replace(trigsimp(Z2), 'Z2', th_2)
    B1 = symo.replace(2*(Z1*Y + Z2*X), 'B1', th_2)
    B2 = symo.replace(2*(Z1*X - Z2*Y), 'B2', th_2)
    B3 = symo.replace(V**2 + W**2 - X**2 - Y**2 - Z1**2 - Z2**2, 'B3', th_2)
    _solve_type_2(symo, [B1, B2, B3], th_2)
    Zi1 = symo.replace(X*cos(th_2) + Y*sin(th_2) + Z1, 'Zi1', th_1)
    Zi2 = symo.replace(X*sin(th_2) - Y*cos(th_2) + Z2, 'Zi2', th_1)
    _solve_type_3(symo, [W, V, Zi1, V, -W, Zi2], th_1)


def _solve_type_8(symo, coeffs, th_1, th_2):
    """Solution for the system:
    X*C1 + Y*C12 = Z1
    X*S1 + Y*S12 = Z2
    """
    X, Y, Z1, Z2 = coeffs
    symo.write_line("# X*cos({0}) + Y*cos({0} + {1}) = Z1".format(th_1, th_2))
    symo.write_line("# X*sin({0}) + Y*sin({0} + {1}) = Z2".format(th_1, th_2))
    X = symo.replace(trigsimp(X), 'X', th_2)
    Y = symo.replace(trigsimp(Y), 'Y', th_2)
    Z1 = symo.replace(trigsimp(Z1), 'Z1', th_2)
    Z2 = symo.replace(trigsimp(Z2), 'Z2', th_2)
    Cj = symo.replace((Z1**2 + Z2**2 - X**2 - Y**2) / (2*X*Y), 'C', th_2)
    YPS = var('YPS%s' % th_2)
    symo.add_to_dict(YPS, (ONE, -ONE))
    symo.add_to_dict(th_2, atan2(YPS*sqrt(1 - Cj**2), Cj))
    Q1 = symo.replace(X + Y*cos(th_2), 'Q1', th_1)
    Q2 = symo.replace(Y*sin(th_2), 'Q2', th_1)
    Den = symo.replace(Q1**2 + Q2**2, 'Den', th_1)
    Si = symo.replace((Q1*Z2 - Q2*Z1)/Den, 'S', th_1)
    Ci = symo.replace((Q1*Z1 + Q2*Z2)/Den, 'C', th_1)
    symo.add_to_dict(th_1, atan2(Si, Ci))


def _solve_square(symo, coeffs, x):
    """ solution for the equation:
    A*x**2 + B*x + C = 0
    """
    A, B, C = coeffs
    A = symo.replace(A, 'A', x)
    B = symo.replace(B, 'B', x)
    C = symo.replace(C, 'C', x)
    Delta = symo.repalce(B**2 - 4*A*C, 'Delta', x)
    YPS = var('YPS%s' % x)
    symo.add_to_dict(YPS, (ONE, -ONE))
    symo.add_to_dict(x, (-B + YPS*sqrt(Delta))/(2*A))








