# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package provides symbolic
modeling of robot dynamics.
"""


from copy import copy

import sympy
from sympy import Matrix

from pysymoro.geometry import compute_screw_transform
from pysymoro.geometry import compute_rot_trans, Transform
from pysymoro.kinematics import compute_vel_acc
from pysymoro.kinematics import compute_omega
from symoroutils import symbolmgr
from symoroutils import tools
from symoroutils.paramsinit import ParamsInit


inert_names = ('XXR', 'XYR', 'XZR', 'YYR', 'YZR',
               'ZZR', 'MXR', 'MYR', 'MZR', 'MR')


def compute_direct_dynamic_NE(robo, symo):
    # antecedent angular velocity, projected into jth frame
    wi = ParamsInit.init_vec(robo)
    w = ParamsInit.init_w(robo)
    jaj = ParamsInit.init_vec(robo, 6)
    # Twist transform list of Matrices 6x6
    jTant = ParamsInit.init_mat(robo, 6)
    beta_star = ParamsInit.init_vec(robo, 6)
    grandJ = ParamsInit.init_mat(robo, 6)
    link_acc = ParamsInit.init_vec(robo, 6)
    H_inv = ParamsInit.init_scalar(robo)
    juj = ParamsInit.init_vec(robo, 6)   # Jj*aj / Hj
    Tau = ParamsInit.init_scalar(robo)
    grandVp = ParamsInit.init_vec(robo, 6)
    grandVp[0] = Matrix([robo.vdot0 - robo.G, robo.w0])
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    for j in xrange(1, robo.NL):
        compute_omega(robo, symo, j, antRj, w, wi)
        compute_screw_transform(robo, symo, j, antRj, antPj, jTant)
        if robo.sigma[j] == 0:
            jaj[j] = Matrix([0, 0, 0, 0, 0, 1])
        elif robo.sigma[j] == 1:
            jaj[j] = Matrix([0, 0, 1, 0, 0, 0])
    for j in xrange(1, robo.NL):
        compute_beta(robo, symo, j, w, beta_star)
        compute_link_acc(robo, symo, j, antRj, antPj, link_acc, w, wi)
        grandJ[j] = inertia_spatial(robo.J[j], robo.MS[j], robo.M[j])
    for j in reversed(xrange(1, robo.NL)):
        replace_beta_J_star(robo, symo, j, grandJ, beta_star)
        compute_Tau(robo, symo, j, grandJ, beta_star, jaj, juj, H_inv, Tau)
        if robo.ant[j] != 0:
            compute_beta_J_star(robo, symo, j, grandJ, jaj, juj, Tau,
                                beta_star, jTant, link_acc)
    for j in xrange(1, robo.NL):
        compute_acceleration(robo, symo, j, jTant, grandVp,
                             juj, H_inv, jaj, Tau, link_acc)
    for j in xrange(1, robo.NL):
        compute_coupled_forces(robo, symo, j, grandVp, grandJ, beta_star)
    return robo.qddot


def compute_beta(robo, symo, j, w, beta_star):
    """
    AVAILABLE: nealgos - compute_beta().
    Internal function. Computes link's wrench when
    the joint accelerations are zero

    Notes
    =====
    beta_star is the output parameter
    """
    E1 = symo.mat_replace(robo.J[j]*w[j], 'JW', j)
    E2 = symo.mat_replace(tools.skew(w[j])*E1, 'KW', j)
    E3 = tools.skew(w[j])*robo.MS[j]
    E4 = symo.mat_replace(tools.skew(w[j])*E3, 'SW', j)
    E5 = -robo.Nex[j] - E2
    E6 = -robo.Fex[j] - E4
    beta_star[j] = Matrix([E6, E5])


def compute_link_acc(robo, symo, j, antRj, antPj, link_acc, w, wi):
    """
    AVAILABLE: nealgos - compute_gamma()
    Internal function. Computes link's accelerations when
    the joint accelerations are zero

    Notes
    =====
    link_acc is the output parameter
    """
    expr1 = tools.skew(wi[j])*Matrix([0, 0, robo.qdot[j]])
    expr1 = symo.mat_replace(expr1, 'WQ', j)
    expr2 = (1 - robo.sigma[j]) * expr1
    expr3 = 2 * robo.sigma[j] * expr1
    expr4 = tools.skew(w[robo.ant[j]]) * antPj[j]
    expr5 = tools.skew(w[robo.ant[j]]) * expr4
    expr6 = antRj[j].transpose() * expr5
    expr7 = expr6 + expr3
    expr7 = symo.mat_replace(expr7, 'LW', j)
    link_acc[j] = Matrix([expr7, expr2])


def replace_beta_J_star(robo, symo, j, grandJ, beta_star):
    """
    AVIALABLE: nealgos - replace_star_terms()
    Internal function. Makes symbol substitution in beta_star
    and grandJ
    """
    grandJ[j] = symo.mat_replace(grandJ[j], 'MJE', j, symmet=True)
    beta_star[j] = symo.mat_replace(beta_star[j], 'VBE', j)


def compute_Tau(robo, symo, j, grandJ, beta_star, jaj, juj, H_inv, Tau):
    """
    SIMILAR: nealgos - compute_tau(), compute_star_terms()
    Internal function. Computes intermediat dynamic variables

    Notes
    =====
    H_inv and Tau are the output parameters
    """
    Jstar_jaj = grandJ[j]*jaj[j]
    if robo.sigma[j] == 2:
        Tau[j] = 0
    else:
        H_inv[j] = 1 / (jaj[j].dot(Jstar_jaj) + robo.IA[j])
        H_inv[j] = symo.replace(H_inv[j], 'JD', j)
        juj[j] = Jstar_jaj*H_inv[j]
        symo.mat_replace(juj[j], 'JU', j)
        joint_friction = robo.fric_s(j) + robo.fric_v(j)
        Tau[j] = jaj[j].dot(beta_star[j]) + robo.GAM[j] - joint_friction
        Tau[j] = symo.replace(Tau[j], 'GW', j)


def compute_beta_J_star(robo, symo, j, grandJ, jaj, juj, Tau,
                        beta_star, jTant, link_acc):
    """
    SIMILAR: nealgos - compute_star_terms()
    Internal function. Computes intermediat dynamic variables

    Notes
    =====
    grandJ and beta_star are the output parameters
    """
    Jstar_jaj = grandJ[j]*jaj[j]
    grandK = symo.mat_replace(grandJ[j] - juj[j]*Jstar_jaj.T,
                              'GK', j)
    E1 = symo.mat_replace(grandK*link_acc[j], 'NG', j)
    E3 = symo.mat_replace(E1 + Tau[j]*juj[j], 'VS', j)
    alpha = symo.mat_replace(E3 - beta_star[j], 'AP', j)
    E4 = symo.mat_replace(jTant[j].T*grandK, 'GX', j)
    E5 = symo.mat_replace(E4*jTant[j], 'TKT', j, symmet=True)
    grandJ[robo.ant[j]] += E5
    beta_star[robo.ant[j]] -= jTant[j].T*alpha


def compute_acceleration(robo, symo, j, jTant, grandVp,
                         juj, H_inv, jaj, Tau, link_acc):
    """
    SIMILAR: nealgos - compute_joint_accel(), compute_link_accel()
    Internal function. Computes joint accelerations and links' twists

    Notes
    =====
    grandVp is the output parameter
    """
    grandR = symo.mat_replace(jTant[j]*grandVp[robo.ant[j]] + link_acc[j],
                              'VR', j)
    E1 = symo.replace(juj[j].dot(grandR), 'GU', j)
    if robo.sigma[j] == 2:
        qddot = 0
    else:
        qddot = H_inv[j]*Tau[j] - E1
    qddot = symo.replace(qddot, 'QDP', j, forced=True)
    grandVp[j] = (grandR + qddot*jaj[j])
    grandVp[j][3:, 0] = symo.mat_replace(grandVp[j][3:, 0], 'WP', j)
    grandVp[j][:3, 0] = symo.mat_replace(grandVp[j][:3, 0], 'VP', j)


def compute_coupled_forces(robo, symo, j, grandVp, grandJ, beta_star):
    """
    AVIALBLE: nealgos - compute_reaction_wrench()
    Internal function.
    """
    E3 = symo.mat_replace(grandJ[j]*grandVp[j], 'DY', j)
    couplforce = E3 - beta_star[j]
    symo.mat_replace(couplforce[3:, 0], 'N', j)
    symo.mat_replace(couplforce[:3, 0], 'E', j)


def inertia_spatial(J, MS, M):
    """
    AVAILABLE: nealgos - inertia_spatial()
    """
    return Matrix([
        (M*sympy.eye(3)).row_join(tools.skew(MS).T),
        tools.skew(MS).row_join(J)
    ])


# TODO:Finish base parameters computation
def base_paremeters(robo_orig):
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
    robo = copy(robo_orig)
    lam = [0 for i in xrange(robo.NL)]
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'regp')
    title = 'Base parameters computation'
    symo.write_params_table(robo, title, inert=True, dynam=True)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    for j in reversed(xrange(1, robo.NL)):
        if robo.sigma[j] == 0:
            # general grouping
            compute_lambda(robo, symo, j, antRj, antPj, lam)
            group_param_rot(robo, symo, j, lam)
            # special grouping
            group_param_rot_spec(robo, symo, j, lam, antRj)
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
        pass
    symo.write_line('*=*')
    symo.write_line()
    title = robo.name + ' grouped inertia parameters'
    symo.write_params_table(robo, title, inert=True, equations=False)
    symo.file_close()
    return robo, symo.sydi


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
    lam03 = lam03.applyfunc(symo.C2S2_simp)
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


def group_param_rot_spec(robo, symo, j, lam, antRj):
    """Internal function. Groups inertia parameters according to the
    special rule for a rotational joint.

    Notes
    =====
    robo is the output paramete
    """
    chainj = robo.chain(j)
    r1, r2, orthog = Transform.find_r12(robo, chainj, antRj, j)
    kRj, all_paral = Transform.kRj(robo, antRj, r1, chainj)
    Kj = robo.get_inert_param(j)
    to_replace = {0, 1, 2, 4, 5, 6, 7}
    if Transform.z_paral(kRj):
        Kj[0] = 0   # XX
        Kj[1] = 0   # XY
        Kj[2] = 0   # XZ
        Kj[4] = 0   # YZ
        to_replace -= {0, 1, 2, 4}
    joint_axis = antRj[chainj[-1]].col(2)
    if all_paral and robo.G.norm() == sympy.Abs(joint_axis.dot(robo.G)):
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
        if Transform.z_paral(kRj):
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
    dotGa = Transform.sna(antRj[j])[2].dot(robo.G)
    if dotGa == tools.ZERO:
        revol_align = robo.ant[robo.ant[j]] == 0 and robo.ant[j] == tools.ZERO
        if robo.ant[j] == 0 or revol_align:
            Kj[9] += robo.IA[j]
            robo.IA[j] = 0
    for i in to_replace:
        Kj[i] = symo.replace(Kj[i], inert_names[i], j)
    robo.put_inert_param(Kj, j)


