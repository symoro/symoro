# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package contains function to compute the base
inertial parameters.
"""


import sympy
from sympy import Matrix

from pysymoro.geometry import compute_rot_trans, Transform
from pysymoro.geometry import Z_AXIS
from symoroutils import tools


inert_names = ('XXR', 'XYR', 'XZR', 'YYR', 'YZR',
               'ZZR', 'MXR', 'MYR', 'MZR', 'MR')


# TODO:Finish base parameters computation
def base_inertial_parameters(robo, symo):
    """Computes grouped inertia parameters. New parametrization
    contains less parameters but generates the same dynamic model

    Parameters
    ==========
    robo : Robot
        Instance of robot description container

    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    lam = [0 for i in xrange(robo.NL)]
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    for j in reversed(xrange(1, robo.NL)):
        if robo.sigma[j] == 0:
            # general grouping
            compute_lambda(robo, symo, j, antRj, antPj, lam)
            group_param_rot(robo, symo, j, lam)
            # special grouping
            group_param_rot_spec(robo, symo, j, lam, antRj, antPj)
            pass
        elif robo.sigma[j] == 1:
            # general grouping
            group_param_prism(robo, symo, j, antRj)
            # special grouping
            group_param_prism_spec(robo, symo, j, antRj, antPj)
            pass
        elif robo.sigma[j] == 2:
            # fixed joint, group everuthing
            compute_lambda(robo, symo, j, antRj, antPj, lam)
            group_param_fix(robo, symo, j, lam)
    symo.write_line('*=*')


def vec_mut_J(v, u):
    """Internal function. Needed for inertia parameters transformation

    Parameters
    ==========
    v, u : Matrix 3x1
        two axis vectors
    Returns : Matrix 6x1
    """
    return Matrix([v[0]*u[0], v[0]*u[1], v[0]*u[2],
                   v[1]*u[1], v[1]*u[2], v[2]*u[2]])


def vec_mut_MS(v, P):
    """Internal function. Needed for inertia parameters transformation

    Parameters
    ==========
    v : Matrix 3x1
        axis vector
    P : Matrix 3x1
        position vector

    Returns : Matrix 6x1
    """
    U = - tools.skew(v)*tools.skew(P)
    return Matrix([2*U[0, 0], U[0, 1] + U[1, 0], U[0, 2] + U[2, 0],
                   2*U[1, 1], U[1, 2] + U[2, 1], 2*U[2, 2]])


def vec_mut_M(P):
    """Internal function. Needed for inertia parameters transformation

    Parameters
    ==========
    P : Matrix 3x1
        position vector

    Returns : Matrix 6x1
    """
    U = -tools.skew(P)*tools.skew(P)
    return Matrix([U[0, 0], U[0, 1], U[0, 2], U[1, 1], U[1, 2], U[2, 2]])


def compute_lambda(robo, symo, j, antRj, antPj, lam):
    """Internal function. Computes the inertia parameters
    transformation matrix

    Notes
    =====
    lam is the output paramete
    """
    lamJJ_list = []
    lamJMS_list = []
    for e1 in xrange(3):
        for e2 in xrange(e1, 3):
            u = vec_mut_J(antRj[j][:, e1], antRj[j][:, e2])
            if e1 != e2:
                u += vec_mut_J(antRj[j][:, e2], antRj[j][:, e1])
            lamJJ_list.append(u.T)
    for e1 in xrange(3):
        v = vec_mut_MS(antRj[j][:, e1], antPj[j])
        lamJMS_list.append(v.T)
    lamJJ = Matrix(lamJJ_list).T  # , 'LamJ', j)
    lamJMS = symo.mat_replace(Matrix(lamJMS_list).T, 'LamMS', j)
    lamJM = symo.mat_replace(vec_mut_M(antPj[j]), 'LamM', j)
    lamJ = lamJJ.row_join(lamJMS).row_join(lamJM)
    lamMS = sympy.zeros(3, 6).row_join(antRj[j]).row_join(antPj[j])
    lamM = sympy.zeros(1, 10)
    lamM[9] = 1
    lam[j] = Matrix([lamJ, lamMS, lamM])


def group_param_rot(robo, symo, j, lam):
    """Internal function. Groups inertia parameters according to the
    general rule for a rotational joint.

    Notes
    =====
    robo is the output paramete
    """
    Kj = robo.get_inert_param(j)

    lam03 = lam[j][:, 0] + lam[j][:, 3]
    lam03 = lam03.applyfunc(tools.C2S2_simp)
    for i in (3, 8, 9):
        Kj[i] = symo.replace(Kj[i], inert_names[i], j)
    if robo.ant[j] != -1:
        Kant = robo.get_inert_param(robo.ant[j])
        Kant += lam03*Kj[3] + lam[j][:, 8]*Kj[8] + lam[j][:, 9]*Kj[9]
        robo.put_inert_param(Kant, robo.ant[j])
    Kj[0] -= Kj[3]  # XX
    Kj[3] = 0   # YY
    Kj[8] = 0   # MZ
    Kj[9] = 0   # M
    robo.put_inert_param(Kj, j)


def group_param_rot_spec(robo, symo, j, lam, antRj, antPj):
    """Internal function. Groups inertia parameters according to the
    special rule for a rotational joint.

    Notes
    =====
    robo is the output paramete
    """
    chainj = robo.chain(j)
    r1, r2, orthog = Transform.find_r12(robo, chainj, antRj, j)
    kRj, all_paral = Transform.kRj(robo, antRj, r1, chainj)
    r1_Px_j, r1_Py_j, r1_Pz_j = Transform.kPj(
        robo, antPj, antRj, r1, chainj
    )
    Kj = robo.get_inert_param(j)
    to_replace = {0, 1, 2, 4, 5, 6, 7}
    if kRj.col(2) == Z_AXIS:
        Kj[0] = 0   # XX
        Kj[1] = 0   # XY
        Kj[2] = 0   # XZ
        Kj[4] = 0   # YZ
        to_replace -= {0, 1, 2, 4}
    joint_axis = antRj[chainj[-1]].col(2)
    if all_paral and \
        (robo.G.norm() == sympy.Abs(joint_axis.dot(robo.G))) and \
        (r1_Px_j == 0) and (r1_Py_j == 0):
        Kj[6] = 0   # MX
        Kj[7] = 0   # MY
        to_replace -= {6, 7}
    if j == r1 or(j == r2 and orthog):
        Kj[5] += robo.IA[j]   # ZZ
        robo.IA[j] = 0
    for i in to_replace:
        Kj[i] = symo.replace(Kj[i], inert_names[i], j)
    robo.put_inert_param(Kj, j)


def group_param_fix(robo, symo, j, lam):
    """Internal function. Groups inertia parameters according to the
    general rule for a fixed joint joint.

    Notes
    =====
    robo is the output paramete
    """
    Kj = robo.get_inert_param(j)
    for i in xrange(10):
        Kj[i] = symo.replace(Kj[i], inert_names[i], j)
    if robo.ant[j] != -1:
        Kant = robo.get_inert_param(robo.ant[j])
        Kant += lam[j]*Kj
        robo.put_inert_param(Kant, robo.ant[j])
    robo.put_inert_param(sympy.zeros(10, 1), j)


def group_param_prism(robo, symo, j, antRj):
    """Internal function. Groups inertia parameters according to the
    general rule for a prismatic joint.

    Notes
    =====
    robo is the output paramete
    """
    Kj = robo.get_inert_param(j)
    for i in xrange(6):
        Kj[i] = symo.replace(Kj[i], inert_names[i], j)
    robo.put_inert_param(Kj, j)
    if robo.ant[j] != -1:
        antJj = antRj[j]*robo.J[j]*antRj[j].T
        robo.J[robo.ant[j]] += antJj
    robo.J[j] = sympy.zeros(3, 3)


def group_param_prism_spec(robo, symo, j, antRj, antPj):
    """Internal function. Groups inertia parameters according to the
    special rule for a prismatic joint.

    Notes
    =====
    robo is the output paramete
    """
    chainj = robo.chain(j)
    r1, r2, orthog = Transform.find_r12(robo, chainj, antRj, j)
    Kj = robo.get_inert_param(j)
    kRj, all_paral = Transform.kRj(robo, antRj, r1, chainj)
    to_replace = {6, 7, 8, 9}
    if r1 < j and j < r2:
        if kRj.col(2) == Z_AXIS:
            Kj[8] = 0   # MZ
            for i in (6, 7):
                Kj[i] = symo.replace(Kj[i], inert_names[i], j)
            robo.MS[robo.ant[j]] += antRj[j]*Matrix([Kj[6], Kj[7], 0])
            robo.JJ[2, 2] -= Kj[6]*antPj[j][0] + Kj[7]*antPj[j][1]
            Kj[6] = 0   # MX
            Kj[7] = 0   # MY
            to_replace -= {6, 7, 8}
        else:
            jar1 = kRj.row(2)
            if jar1[2] != 0:
                Kj[6] -= jar1[0]/jar1[2]*Kj[8]
                Kj[7] -= jar1[1]/jar1[2]*Kj[8]
                Kj[8] = 0   # MZ
                to_replace -= {8}
            elif jar1[0]*jar1[1] != 0:
                Kj[6] -= jar1[0]/jar1[1]*Kj[7]
                Kj[7] = 0   # MY
                to_replace -= {7}
            elif jar1[0] != 0:
                Kj[7] = 0   # MY
                to_replace -= {7}
            else:
                Kj[6] = 0   # MX
                to_replace -= {6}
    elif j < r1:
        Kj[6] = 0   # MX
        Kj[7] = 0   # MY
        Kj[8] = 0   # MZ
        to_replace -= {6, 7, 8}
    #TOD: rewrite
    dotGa = antRj[j].col(2).dot(robo.G)
    if dotGa == tools.ZERO:
        revol_align = robo.ant[robo.ant[j]] == 0 and robo.ant[j] == tools.ZERO
        if robo.ant[j] == 0 or revol_align:
            Kj[9] += robo.IA[j]
            robo.IA[j] = 0
    for i in to_replace:
        Kj[i] = symo.replace(Kj[i], inert_names[i], j)
    robo.put_inert_param(Kj, j)


