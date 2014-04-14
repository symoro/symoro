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
    return Matrix([[0, -vec[2], vec[1]],
                   [vec[2], 0, -vec[0]],
                   [-vec[1], vec[0], 0]])


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


def get_max_coef_mul(sym, x_term):
    k, ex = x_term.as_coeff_Mul()
    coef = sym / k
    pow_x = ex.as_powers_dict()
    pow_c = coef.as_powers_dict()
    pow_c[-1] = 0
    for j, pow_j in pow_x.iteritems():
        num_j = -j
        if j in pow_c and pow_c[j] >= pow_j:
            pow_c[j] -= pow_j
        elif num_j in pow_c and pow_c[num_j] >= pow_j:
            pow_c[num_j] -= pow_j
            if pow_j % 2:
                pow_c[-1] += 1
        else:
            return ZERO
    return Mul.fromiter(c**p for c, p in pow_c.iteritems())


def get_max_coef_list(sym, x_term):
    return [get_max_coef_mul(s, x_term) for s in Add.make_args(sym)]


def get_max_coef(sym, x_term):
    return Add.fromiter(
        get_max_coef_mul(s, x_term) for s in Add.make_args(sym)
    )


def get_pos_neg(str_term):
    if str_term.find('m') != -1:
        s_split = str_term.split('m')
        return s_split[0], s_split[1]
    else:
        return str_term, ''


def reduce_str(str1, str2):
    while True:
        for j, char in enumerate(str1):
            if char in 'AG':
                i = str2.find(str1[j:j+2])
                k = 2
            else:
                i = str2.find(char)
                k = 1
            if i != -1:
                if i+k < len(str2):
                    str2_tail = str2[i+k:]
                else:
                    str2_tail = ''
                if j+k < len(str1):
                    str1_tail = str1[j+k:]
                else:
                    str1_tail = ''
                str2 = str2[:i] + str2_tail
                str1 = str1[:j] + str1_tail
                break
        else:
            break
    return str1, str2


def ang_sum(np1, np2, nm1, nm2):
    np2, nm1 = reduce_str(np2, nm1)
    np1, nm2 = reduce_str(np1, nm2)
    if len(nm1) + len(nm2) == 0:
        return np1 + np2
    else:
        return np1 + np2 + 'm' + nm1 + nm2


def cos_sin_syms(name):
    if isinstance(name, str) and name[0] == 'm':
        cos_term, sin_term = var('C{0}, S{0}'.format(name[1:]))
        return cos_term, -sin_term
    else:
        return var('C{0}, S{0}'.format(name))


def sym_less(val_a, val_b):
    val_a_measure = val_a.count_ops()
    val_b_measure = val_b.count_ops()
    return val_a_measure < val_b_measure


def get_angles(expr):
    angles_s = set()
    for sin_term in expr.atoms(sin):
        angles_s |= set(sin_term.args)
    angles_c = set()
    for cos_term in expr.atoms(cos):
        angles_c |= set(cos_term.args)
    return angles_s & angles_c


def cancel_terms(sym, x_term, coef):
    if coef.is_Add:
        for arg_c in coef.args:
            sym = cancel_terms(sym, x_term, arg_c)
    else:
        terms = Add.make_args(sym)
        return Add.fromiter(t for t in terms if t != x_term*coef)


def trignometric_info(sym):
    if not sym.has(sin) and not sym.has(cos):
        short_form = True
        names = get_trig_couple_names(sym)
    else:
        short_form = False
        names = get_angles(sym)
    return names, short_form


