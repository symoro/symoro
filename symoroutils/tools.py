# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module contains the miscellaneous helper functions needed by the
SYMORO software package.
"""

import itertools
import re

from sympy import Matrix, Expr, Integer
from sympy import sin, cos
from sympy import Mul, Add


ZERO = Integer(0)
ONE = Integer(1)
FAIL = 1
OK = 0
CLOSED_LOOP = 'Closed loop'
SIMPLE = 'Simple'
TREE = 'Tree'
TYPES = [SIMPLE, TREE, CLOSED_LOOP]
INT_KEYS = ['ant', 'sigma', 'mu']


def skew(vec):
    """skew-symmetry : Generates vectorial preproduct matrix

    Parameters
    ==========
    v: Matrix 3x1
        vector

    Returns
    =======
    hat: Matrix 3x3
    """
    return Matrix([
        [0, -vec[2], vec[1]],
        [vec[2], 0, -vec[0]],
        [-vec[1], vec[0], 0]
    ])


def skew2vec(mat):
    """
    Return a 3x1 vector from 3x3 skew-symmetric matrix.

    Args:
        mat: A 3x3 skew-symmetric matrix (Matrix)
    Returns:
        A 3x1 vector (Matrix)
    """
    vec0 = mat[2, 1]
    vec1 = mat[0, 2]
    vec2 = mat[1, 0]
    return Matrix([vec0, vec1, vec2])


def l2str(list_var, spacing=8):
    """Converts a list into string, that will be
    written into the text table.

    Parameters
    ==========
    list_var: list
        List to be converted
    spacing: int, optional
        Defines the size of one cell of the table

    Returns
    =======
    ret_str: string
        String representation

    Notes
    =====
    l2str([1, 2, 3]) will be converted into '1      2      3      '
    """
    ret_str = ''
    for i in list_var:
        ret_str += str(i) + ' '*(spacing-len(str(i)))
    return ret_str


def find_trig_names(sym, pref=r'', pref_len=0, post=r'', post_len=0):
    search_res = re.findall(pref + r'[AGm0-9]*' + post, str(sym))
    if post_len == 0:
        return set([s[pref_len:] for s in search_res])
    else:
        return set([s[pref_len:-post_len] for s in search_res])


def get_trig_couple_names(sym):
    names_s = find_trig_names(sym, r'S', 1)
    names_c = find_trig_names(sym, r'C', 1)
    return names_c & names_s


#TODO get rid of it
def get_max_coef_mul(sym, x_term):
    pow_x = x_term.as_powers_dict()
    pow_c = sym.as_powers_dict()
    for x, pow_j in pow_x.iteritems():
        if x in pow_c and pow_c[x] >= pow_j:
            pow_c[x] -= pow_j
        else:
            return None
    return Mul.fromiter(c**p for c, p in pow_c.iteritems())


def factorize(sym, x_term):
    """sym = rest + coeff*x_term"""
    rest = ZERO
    coeff = ZERO
    for s in Add.make_args(sym):
        k = get_max_coef_mul(s, x_term)
        if k is None:
            rest += s
        else:
            coeff += k
    return rest, coeff


def get_coeffs(sym, var_list):
    coeffs = []
    for x in var_list:
        sym, c = factorize(sym, x)
        print
        print x
        print sym
        print c
        coeffs.append(c)
    coeffs.append(sym)
    return coeffs


def get_angles(expr):
    angles_s = set()
    for sin_term in expr.atoms(sin):
        angles_s |= set(sin_term.args)
    angles_c = set()
    for cos_term in expr.atoms(cos):
        angles_c |= set(cos_term.args)
    return angles_s & angles_c


def simplify(sym, C2S2=False):
    if isinstance(sym, Matrix):
        for i in xrange(len(sym)):
            sym[i] = simplify(sym[i])
        return sym
    elif isinstance(sym, Expr):
        if isinstance(sym, Add):
            sym = Add.fromiter(simplify(x) for x in sym.args)
            if C2S2:
                sym = C2S2_simp(sym)
            sym = CS12_simp(sym)
            return sym
        else:
            new_syms = dict()
            for arg in sym.args:
                if arg.has(Add) and arg.has(sin) and arg.has(cos):
                    new_syms[arg] = simplify(arg)
            return sym.xreplace(new_syms)


def C2S2_simp(sym):
    """
    Example
    =======
    >> print C2S2_simp(sympify("-C**2*RL + S*(D - RL*S)"))
    D*S - RL
    """
    if not sym.is_Add:
        repl_dict = {}
        for term in sym.atoms(Add):
            repl_dict[term] = C2S2_simp(term)
        sym = sym.xreplace(repl_dict)
        return sym
    else:
        sym_old = sym
        while True:
            angles = get_angles(sym)
            for x in angles:
                if not sym.is_Add:
                    break
                sym = try_opt(ONE, None, cos(x)**2, sin(x)**2, sym)
            if sym_old == sym:
                break
            else:
                sym_old = sym
        return sym


def CS12_simp(sym):
    """
    Example
    =======
    >> print SymbolManager().CS12_simp(sympify("C2*C3 - S2*S3"))
    C23 = C2*C3 - S2*S3
    C23
    >> print SymbolManager().CS12_simp(sympify("C2*S3*R + S2*C3*R"))
    S23 = C2*S3 + S2*C3
    R*S23
    """
    if not sym.is_Add:
        repl_dict = {}
        for term in sym.atoms(Add):
            repl_dict[term] = CS12_simp(term)
        sym = sym.xreplace(repl_dict)
        return sym
    else:
        names = get_angles(sym)
        res = sym
        for n1, n2 in itertools.combinations(names, 2):
            if not res.is_Add:
                break
            C1, S1 = cos(n1), sin(n1)
            C2, S2 = cos(n2), sin(n2)
            C12, S12 = cos(n1+n2), sin(n1+n2)
            C1m2, S1m2 = cos(n1-n2), sin(n1-n2)
            res = try_opt(S12, S1m2, S1*C2, C1*S2, res)
            if not res.is_Add:
                break
            res = try_opt(C1m2, C12, C1*C2, S1*S2, res)
        if res != sym:
            return CS12_simp(res)
        else:
            return sym


def is_less(A, B):
    return A.count_ops() < B.count_ops


def try_opt(A, Am, B, C, old_sym, silent=False):
    """Replaces B + C by A or B - C by Am.
    Chooses the best option.
    """
    assert isinstance(old_sym, Add)

    _, KB = factorize(old_sym, B)
    _, KC = factorize(old_sym, C)
    if KB == ZERO or KC == ZERO:
        return old_sym
    else:
        KB_args = Add.make_args(KB)
        KC_args = set(Add.make_args(KC))
        res = old_sym
        for k in KB_args:
            if k in KC_args:
                res = res - k*B - k*C + k*A
            elif -k in KC_args:
                res = res - k*B + k*C + k*Am
        return res


#from sympy import sympify, Wild, sin, cos, factor
#
#eq = sympify("""D3**2*RL4*sin(th1)**6*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th2) + D3**2*RL4*sin(th1)**6*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th2)*cos(th2 + th3)**2 + D3**2*RL4*sin(th1)**6*sin(th2)*sin(th5)*sin(th2 + th3)**3*cos(th2)*cos(th4)**2 + D3**2*RL4*sin(th1)**6*sin(th2)*sin(th5)*sin(th2 + th3)*cos(th2)*cos(th4)**2*cos(th2 + th3)**2 + D3**2*RL4*sin(th1)**6*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th2)**2*cos(th2 + th3) + D3**2*RL4*sin(th1)**6*sin(th4)**2*sin(th5)*cos(th2)**2*cos(th2 + th3)**3 + D3**2*RL4*sin(th1)**6*sin(th5)*sin(th2 + th3)**2*cos(th2)**2*cos(th4)**2*cos(th2 + th3) +
#D3**2*RL4*sin(th1)**6*sin(th5)*cos(th2)**2*cos(th4)**2*cos(th2 + th3)**3 + 3*D3**2*RL4*sin(th1)**4*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**2*cos(th2) + 3*D3**2*RL4*sin(th1)**4*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th1)**2*cos(th2)*cos(th2 + th3)**2 + 3*D3**2*RL4*sin(th1)**4*sin(th2)*sin(th5)*sin(th2 + th3)**3*cos(th1)**2*cos(th2)*cos(th4)**2 + 3*D3**2*RL4*sin(th1)**4*sin(th2)*sin(th5)*sin(th2 + th3)*cos(th1)**2*cos(th2)*cos(th4)**2*cos(th2 + th3)**2 + 3*D3**2*RL4*sin(th1)**4*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**2*cos(th2)**2*cos(th2 + th3) +
#3*D3**2*RL4*sin(th1)**4*sin(th4)**2*sin(th5)*cos(th1)**2*cos(th2)**2*cos(th2 + th3)**3 + 3*D3**2*RL4*sin(th1)**4*sin(th5)*sin(th2 + th3)**2*cos(th1)**2*cos(th2)**2*cos(th4)**2*cos(th2 + th3) + 3*D3**2*RL4*sin(th1)**4*sin(th5)*cos(th1)**2*cos(th2)**2*cos(th4)**2*cos(th2 + th3)**3 + 3*D3**2*RL4*sin(th1)**2*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**4*cos(th2) + 3*D3**2*RL4*sin(th1)**2*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th1)**4*cos(th2)*cos(th2 + th3)**2 +
#3*D3**2*RL4*sin(th1)**2*sin(th2)*sin(th5)*sin(th2 + th3)**3*cos(th1)**4*cos(th2)*cos(th4)**2 + 3*D3**2*RL4*sin(th1)**2*sin(th2)*sin(th5)*sin(th2 + th3)*cos(th1)**4*cos(th2)*cos(th4)**2*cos(th2 + th3)**2 + 3*D3**2*RL4*sin(th1)**2*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**4*cos(th2)**2*cos(th2 + th3) + 3*D3**2*RL4*sin(th1)**2*sin(th4)**2*sin(th5)*cos(th1)**4*cos(th2)**2*cos(th2 + th3)**3 + 3*D3**2*RL4*sin(th1)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**4*cos(th2)**2*cos(th4)**2*cos(th2 + th3) + 3*D3**2*RL4*sin(th1)**2*sin(th5)*cos(th1)**4*cos(th2)**2*cos(th4)**2*cos(th2 + th3)**3 +
#D3**2*RL4*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**6*cos(th2) + D3**2*RL4*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th1)**6*cos(th2)*cos(th2 + th3)**2 + D3**2*RL4*sin(th2)*sin(th5)*sin(th2 + th3)**3*cos(th1)**6*cos(th2)*cos(th4)**2 + D3**2*RL4*sin(th2)*sin(th5)*sin(th2 + th3)*cos(th1)**6*cos(th2)*cos(th4)**2*cos(th2 + th3)**2 + D3**2*RL4*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**6*cos(th2)**2*cos(th2 + th3) + D3**2*RL4*sin(th4)**2*sin(th5)*cos(th1)**6*cos(th2)**2*cos(th2 + th3)**3 + D3**2*RL4*sin(th5)*sin(th2 + th3)**2*cos(th1)**6*cos(th2)**2*cos(th4)**2*cos(th2 + th3) +
#D3**2*RL4*sin(th5)*cos(th1)**6*cos(th2)**2*cos(th4)**2*cos(th2 + th3)**3 - D3*RL4**2*sin(th1)**6*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**4 - D3*RL4**2*sin(th1)**6*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th2 + th3)**2 - D3*RL4**2*sin(th1)**6*sin(th2)*sin(th5)*sin(th2 + th3)**4*cos(th4)**2 - D3*RL4**2*sin(th1)**6*sin(th2)*sin(th5)*sin(th2 + th3)**2*cos(th4)**2*cos(th2 + th3)**2 - D3*RL4**2*sin(th1)**6*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th2)*cos(th2 + th3) - D3*RL4**2*sin(th1)**6*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th2)*cos(th2 + th3)**3 - D3*RL4**2*sin(th1)**6*sin(th5)*sin(th2 + th3)**3*cos(th2)*cos(th4)**2*cos(th2 + th3) -
#D3*RL4**2*sin(th1)**6*sin(th5)*sin(th2 + th3)*cos(th2)*cos(th4)**2*cos(th2 + th3)**3 - 3*D3*RL4**2*sin(th1)**4*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**4*cos(th1)**2 - 3*D3*RL4**2*sin(th1)**4*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**2*cos(th2 + th3)**2 - 3*D3*RL4**2*sin(th1)**4*sin(th2)*sin(th5)*sin(th2 + th3)**4*cos(th1)**2*cos(th4)**2 - 3*D3*RL4**2*sin(th1)**4*sin(th2)*sin(th5)*sin(th2 + th3)**2*cos(th1)**2*cos(th4)**2*cos(th2 + th3)**2 - 3*D3*RL4**2*sin(th1)**4*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**2*cos(th2)*cos(th2 + th3) -
#3*D3*RL4**2*sin(th1)**4*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th1)**2*cos(th2)*cos(th2 + th3)**3 - 3*D3*RL4**2*sin(th1)**4*sin(th5)*sin(th2 + th3)**3*cos(th1)**2*cos(th2)*cos(th4)**2*cos(th2 + th3) - 3*D3*RL4**2*sin(th1)**4*sin(th5)*sin(th2 + th3)*cos(th1)**2*cos(th2)*cos(th4)**2*cos(th2 + th3)**3 - 3*D3*RL4**2*sin(th1)**2*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**4*cos(th1)**4 - 3*D3*RL4**2*sin(th1)**2*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**4*cos(th2 + th3)**2 - 3*D3*RL4**2*sin(th1)**2*sin(th2)*sin(th5)*sin(th2 + th3)**4*cos(th1)**4*cos(th4)**2 -
#3*D3*RL4**2*sin(th1)**2*sin(th2)*sin(th5)*sin(th2 + th3)**2*cos(th1)**4*cos(th4)**2*cos(th2 + th3)**2 - 3*D3*RL4**2*sin(th1)**2*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**4*cos(th2)*cos(th2 + th3) - 3*D3*RL4**2*sin(th1)**2*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th1)**4*cos(th2)*cos(th2 + th3)**3 - 3*D3*RL4**2*sin(th1)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**4*cos(th2)*cos(th4)**2*cos(th2 + th3) - 3*D3*RL4**2*sin(th1)**2*sin(th5)*sin(th2 + th3)*cos(th1)**4*cos(th2)*cos(th4)**2*cos(th2 + th3)**3 - D3*RL4**2*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**4*cos(th1)**6 - D3*RL4**2*sin(th2)*sin(th4)**2*sin(th5)*sin(th2 + th3)**2*cos(th1)**6*cos(th2 + th3)**2 -
#D3*RL4**2*sin(th2)*sin(th5)*sin(th2 + th3)**4*cos(th1)**6*cos(th4)**2 - D3*RL4**2*sin(th2)*sin(th5)*sin(th2 + th3)**2*cos(th1)**6*cos(th4)**2*cos(th2 + th3)**2 - D3*RL4**2*sin(th4)**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**6*cos(th2)*cos(th2 + th3) - D3*RL4**2*sin(th4)**2*sin(th5)*sin(th2 + th3)*cos(th1)**6*cos(th2)*cos(th2 + th3)**3 - D3*RL4**2*sin(th5)*sin(th2 + th3)**3*cos(th1)**6*cos(th2)*cos(th4)**2*cos(th2 + th3) - D3*RL4**2*sin(th5)*sin(th2 + th3)*cos(th1)**6*cos(th2)*cos(th4)**2*cos(th2 + th3)**3""")
#
#print eq
#
#eq = simplify(eq, C2S2=True)
#print eq
#print factor(eq)
#print "###"

#print get_coeffs(sympify('X*sin(theta1) - Y + Z*sin(theta1)'),
#           [sympify('sin(theta1)'), sympify('cos(theta1)')])
