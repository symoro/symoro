# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module contains the functions for the computation of Inertia
matrix.
"""


import sympy
from sympy import Matrix

from pysymoro.geometry import compute_rot_trans
from symoroutils.paramsinit import ParamsInit
from symoroutils import tools


CHARSYMS = (
    'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'j', 'k', 'l', 'm', 'n',
    'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', 'aa', 'ab',
    'ac', 'ad', 'ae', 'af', 'ag', 'ah', 'aj', 'ak', 'al', 'am', 'an',
    'ap', 'aq', 'ar', 'as', 'at', 'au', 'av', 'aw', 'ax', 'ay', 'az',
    'ba', 'bb', 'bc', 'bd', 'be', 'bf', 'bg', 'bh', 'bj', 'bk', 'bl',
    'bm', 'bn', 'bp', 'bq', 'br', 'bs', 'bt', 'bu', 'bv', 'bw', 'bx',
    'by', 'bz'
)


def inertia_spatial(inertia, ms_tensor, mass):
    """
    Setup spatial inertia matrix (internal function).
    """
    return Matrix([
        (mass * sympy.eye(3)).row_join(tools.skew(ms_tensor).transpose()),
        tools.skew(ms_tensor).row_join(inertia)
    ])


def replace_composite_terms(
    symo, j, comp_inertia3, comp_ms, comp_mass
):
    """
    Replace composite inertia terms (internal function).

    Note:
        comp_inertia3, comp_ms, comp_mass are the output parameters.
    """
    comp_inertia3[j] = symo.mat_replace(comp_inertia3[j], 'JP', j)
    comp_ms[j] = symo.mat_replace(comp_ms[j], 'MSP', j)
    comp_mass[j] = symo.replace(comp_mass[j], 'MP', j)


def compute_composite_inertia(
    robo, symo, j, antRj, antPj,
    aje1, comp_inertia3, comp_ms, comp_mass
):
    """
    Compute composite inertia (internal function).

    Note:
        aje1, comp_inertia3, comp_ms, comp_mass are the output
        parameters.
    """
    i = robo.ant[j]
    i_ms_j_c = antRj[j] * comp_ms[j]
    i_ms_j_c = symo.mat_replace(i_ms_j_c, 'AS', j)
    expr1 = antRj[j] * comp_inertia3[j]
    expr1 = symo.mat_replace(expr1, 'AJ', j)
    aje1[j] = expr1[:, 2]
    expr2 = expr1 * antRj[j].transpose()
    expr2 = symo.mat_replace(expr2, 'AJA', j)
    expr3 = tools.skew(antPj[j]) * tools.skew(i_ms_j_c)
    expr3 = symo.mat_replace(expr3, 'PAS', j)
    comp_inertia3[i] += expr2 - (expr3 + expr3.transpose()) + \
        (comp_mass[j] * tools.skew(antPj[j]) * \
        tools.skew(antPj[j]).transpose())
    comp_ms[i] = comp_ms[i] + i_ms_j_c + (antPj[j] * comp_mass[j])
    comp_mass[i] = comp_mass[i] + comp_mass[j]


def compute_diagonal_elements(
    robo, symo, j, comp_inertia3, comp_ms, comp_mass,
    forces, moments, inertia_a22
):
    """
    Compute diagonal elements of the inertia matrix (internal function).

    Note:
        forces, moments, inertia_a22 are the output parameters
    """
    if robo.sigma[j] == 0:
        forces[j] = Matrix([-comp_ms[j][1], comp_ms[j][0], 0])
        moments[j] = comp_inertia3[j][:, 2]
        inertia_a22[j-1, j-1] = comp_inertia3[j][2, 2] + robo.IA[j]
    elif robo.sigma[j] == 1:
        forces[j] = Matrix([0, 0, comp_mass[j]])
        moments[j] = Matrix([comp_ms[j][1], -comp_ms[j][0], 0])
        inertia_a22[j-1, j-1] = comp_mass[j] + robo.IA[j]
    forces[j] = symo.mat_replace(forces[j], 'E' + CHARSYMS[j], j)
    moments[j] = symo.mat_replace(moments[j], 'N' + CHARSYMS[j], j)


def compute_triangle_elements(
    robo, symo, j, k, ka, antRj, antPj, aje1,
    forces, moments, inertia_a12, inertia_a22
):
    """
    Compute elements below and above diagonal of the inertia matrix
    (internal function).

    Note:
        forces, moments, inertia_a12, inertia_a22 are the output
        parameters
    """
    forces[ka] = antRj[k] * forces[k]
    if k == j and robo.sigma[j] == 0:
        moments[ka] = aje1[k] + \
            (tools.skew(antPj[k]) * antRj[k] * forces[k])
    else:
        moments[ka] = (antRj[k] * moments[k]) + \
            (tools.skew(antPj[k]) * antRj[k] * forces[k])
    if ka == 0:
        inertia_a12[j][:3, 0] = symo.mat_replace(
            forces[ka], 'AV0', j, forced=True
        )
        inertia_a12[j][3:, 0] = symo.mat_replace(
            moments[ka], 'AW0', j, forced=True
        )
    else:
        symo.mat_replace(forces[ka], 'E' + CHARSYMS[j], ka)
        symo.mat_replace(moments[ka], 'N' + CHARSYMS[j], ka)
        if robo.sigma[ka] == 0:
            inertia_a22[j-1, ka-1] = moments[ka][2]
        elif robo.sigma[ka] == 1:
            inertia_a22[j-1, ka-1] = forces[ka][2]
        inertia_a22[ka-1, j-1] = inertia_a22[j-1, ka-1]


def fixed_inertia_matrix(robo, symo):
    """
    Compute Inertia Matrix for robots with fixed base. This function
    computes just the A22 matrix when the inertia matrix
    A = [A11, A12; A12.transpose(), A22].
    """
    # init terms
    comp_inertia3, comp_ms, comp_mass = ParamsInit.init_jplus(robo)
    aje1 = ParamsInit.init_vec(robo)
    forces = ParamsInit.init_vec(robo, ext=1)
    moments = ParamsInit.init_vec(robo, ext=1)
    inertia_a12 = ParamsInit.init_vec(robo, num=6)
    inertia_a22 = sympy.zeros(robo.nl, robo.nl)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    for j in reversed(range(1, robo.NL)):
        replace_composite_terms(
            symo, j, comp_inertia3, comp_ms, comp_mass
        )
        if j != 1:
            compute_composite_inertia(
                robo, symo, j, antRj, antPj,
                aje1, comp_inertia3, comp_ms, comp_mass
            )
    for j in range(1, robo.NL):
        compute_diagonal_elements(
            robo, symo, j, comp_inertia3, comp_ms,
            comp_mass, forces, moments, inertia_a22
        )
        ka = j
        while ka != 1:
            k = ka
            ka = robo.ant[ka]
            compute_triangle_elements(
                robo, symo, j, k, ka, antRj, antPj, aje1,
                forces, moments, inertia_a12, inertia_a22
            )
    symo.mat_replace(inertia_a22, 'A', forced=True, symmet=True)
    return inertia_a22


def floating_inertia_matrix(robo, symo):
    """
    Compute Inertia Matrix for robots with floating or mobile base. This
    function computes the A11, A12 and A22 matrices when the inertia
    matrix A = [A11, A12; A12.transpose(), A22]
    """
    # init terms
    comp_inertia3, comp_ms, comp_mass = ParamsInit.init_jplus(robo)
    aje1 = ParamsInit.init_vec(robo)
    forces = ParamsInit.init_vec(robo, ext=1)
    moments = ParamsInit.init_vec(robo, ext=1)
    inertia_a12 = ParamsInit.init_vec(robo, num=6)
    inertia_a22 = sympy.zeros(robo.nl, robo.nl)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    for j in reversed(range(0, robo.NL)):
        replace_composite_terms(
            symo, j, comp_inertia3, comp_ms, comp_mass
        )
        if j != 0:
            compute_composite_inertia(
                robo, symo, j, antRj, antPj,
                aje1, comp_inertia3, comp_ms, comp_mass
            )
    for j in range(1, robo.NL):
        compute_diagonal_elements(
            robo, symo, j, comp_inertia3, comp_ms,
            comp_mass, forces, moments, inertia_a22
        )
        ka = j
        while ka != 0:
            k = ka
            ka = robo.ant[ka]
            compute_triangle_elements(
                robo, symo, j, k, ka, antRj, antPj, aje1,
                forces, moments, inertia_a12, inertia_a22
            )
    symo.mat_replace(inertia_a22, 'A', forced=True, symmet=True)
    inertia_a11 = inertia_spatial(
        comp_inertia3[0], comp_ms[0], comp_mass[0]
    )
    inertia_a11 = symo.mat_replace(
        inertia_a11, 'Jcomp', 0, forced=True, symmet=True
    )
    # setup inertia_a12 in Matrix form
    a12mat = sympy.zeros(6, robo.NL)
    for j in range(1, robo.NL):
        a12mat[:, j] = inertia_a12[j]
    a12mat = a12mat[:, 1:]
    # setup the complete inertia matrix
    inertia = Matrix([
        inertia_a11.row_join(a12mat),
        a12mat.transpose().row_join(inertia_a22)
    ])
    return inertia


