# -*- coding: utf-8 -*-


"""
This module contains the miscellaneous helper functions needed by the
SYMORO software package.
"""


import re

from sympy import Matrix
from sympy import Integer
from sympy import sin, cos
from sympy import Mul, Add, var


ZERO = Integer(0)
ONE = Integer(1)
FAIL = 1
OK = 0
CLOSED_LOOP = 'Closed loop'
SIMPLE = 'Simple'
TREE = 'Tree'
TYPES = [SIMPLE, TREE, CLOSED_LOOP]
INT_KEYS = ['ant', 'sigma', 'mu']


def skew(v):
    """skew-symmetry : Generates vectorial preproduct matrix

    Parameters
    ==========
    v: Matrix 3x1
        vector

    Returns
    =======
    hat: Matrix 3x3
    """
    return Matrix([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])


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
    s: string
        String representation

    Notes
    =====
    l2str([1, 2, 3]) will be converted into '1      2      3      '
    """
    s = ''
    for i in list_var:
        s += str(i) + ' '*(spacing-len(str(i)))
    return s


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


def get_max_coef_mul(sym, x):
    k, ex = x.as_coeff_Mul()
    coef = sym / k
    pow_x = ex.as_powers_dict()
    pow_c = coef.as_powers_dict()
    pow_c[-1] = 0
    for a, pa in pow_x.iteritems():
        na = -a
        if a in pow_c and pow_c[a] >= pa:
            pow_c[a] -= pa
        elif na in pow_c and pow_c[na] >= pa:
            pow_c[na] -= pa
            if pa % 2:
                pow_c[-1] += 1
        else:
            return ZERO
    return Mul.fromiter(c**p for c, p in pow_c.iteritems())


def get_max_coef_list(sym, x):
    return [get_max_coef_mul(s, x) for s in Add.make_args(sym)]


def get_max_coef(sym, x):
    return Add.fromiter(get_max_coef_mul(s, x) for s in Add.make_args(sym))


def get_pos_neg(s):
    if s.find('m') != -1:
        s_split = s.split('m')
        return s_split[0], s_split[1]
    else:
        return s, ''


def reduce_str(s1, s2):
    while True:
        for j, char in enumerate(s1):
            if char in 'AG':
                i = s2.find(s1[j:j+2])
                k = 2
            else:
                i = s2.find(char)
                k = 1
            if i != -1:
                if i+k < len(s2):
                    s2_tail = s2[i+k:]
                else:
                    s2_tail = ''
                if j+k < len(s1):
                    s1_tail = s1[j+k:]
                else:
                    s1_tail = ''
                s2 = s2[:i] + s2_tail
                s1 = s1[:j] + s1_tail
                break
        else:
            break
    return s1, s2


def ang_sum(np1, np2, nm1, nm2):
    np2, nm1 = reduce_str(np2, nm1)
    np1, nm2 = reduce_str(np1, nm2)
    if len(nm1) + len(nm2) == 0:
        return np1 + np2
    else:
        return np1 + np2 + 'm' + nm1 + nm2


def CS_syms(name):
    if isinstance(name, str) and name[0] == 'm':
        C, S = var('C{0}, S{0}'.format(name[1:]))
        return C, -S
    else:
        return var('C{0}, S{0}'.format(name))


def sym_less(A, B):
    A_measure = A.count_ops()
    B_measure = B.count_ops()
    return A_measure < B_measure


def get_angles(expr):
    angles_s = set()
    for s in expr.atoms(sin):
        angles_s |= set(s.args)
    angles_c = set()
    for c in expr.atoms(cos):
        angles_c |= set(c.args)
    return angles_s & angles_c


def cancel_terms(sym, X, coef):
    if coef.is_Add:
        for arg_c in coef.args:
            sym = cancel_terms(sym, X, arg_c)
    else:
        terms = Add.make_args(sym)
        return Add.fromiter(t for t in terms if t != X*coef)


def trignometric_info(sym):
    if not sym.has(sin) and not sym.has(cos):
        short_form = True
        names = get_trig_couple_names(sym)
    else:
        short_form = False
        names = get_angles(sym)
    return names, short_form


