#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module of SYMORO package provides symbolic
modeling of robot dynamics.
"""


import sympy
from sympy import Matrix
from copy import copy, deepcopy

from pysymoro.geometry import compute_screw_transform
from pysymoro.geometry import compute_rot_trans, Transform
from pysymoro.kinematics import compute_vel_acc
from pysymoro.kinematics import compute_omega
from symoroutils import symbolmgr
from symoroutils import tools
from symoroutils.paramsinit import ParamsInit


chars = 'ABCDEFGHJKLMNPQRSTUVWXYZ'
inert_names = ('XXR', 'XYR', 'XZR', 'YYR', 'YZR',
               'ZZR', 'MXR', 'MYR', 'MZR', 'MR')


def Newton_Euler(robo, symo):
    """Internal function. Computes Inverse Dynamic Model using
    Newton-Euler formulation

    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    symo : symbolmgr.SymbolManager
        Instance of symbolic manager
    """
    # init external forces
    Fex = copy(robo.Fex)
    Nex = copy(robo.Nex)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    # init velocities and accelerations
    w, wdot, vdot, U = compute_vel_acc(robo, symo, antRj, antPj)
    # init forces vectors
    F = ParamsInit.init_vec(robo)
    N = ParamsInit.init_vec(robo)
    Fjnt = ParamsInit.init_vec(robo)
    Njnt = ParamsInit.init_vec(robo)
    for j in xrange(1, robo.NL):
        compute_wrench(robo, symo, j, w, wdot, U, vdot, F, N)
    for j in reversed(xrange(1, robo.NL)):
        compute_joint_wrench(robo, symo, j, antRj, antPj, vdot,
                             Fjnt, Njnt, F, N, Fex, Nex)
    for j in xrange(1, robo.NL):
        compute_torque(robo, symo, j, Fjnt, Njnt)


def dynamic_identification_NE(robo):
    """Computes Dynamic Identification Model using
    Newton-Euler formulation

    Parameters
    ==========
    robo : Robot
        Instance of robot description container

    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """

    # init forces vectors
    Fjnt = ParamsInit.init_vec(robo)
    Njnt = ParamsInit.init_vec(robo)
    # init file output, writing the robot description
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'dim')
    title = "Dynamic identification model using Newton - Euler Algorith"
    symo.write_params_table(robo, title, inert=True, dynam=True)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    # init velocities and accelerations
    w, wdot, vdot, U = compute_vel_acc(robo, symo, antRj, antPj)
    # virtual robot with only one non-zero parameter at once
    robo_tmp = deepcopy(robo)
    robo_tmp.IA = sympy.zeros(robo.NL, 1)
    robo_tmp.FV = sympy.zeros(robo.NL, 1)
    robo_tmp.FS = sympy.zeros(robo.NL, 1)
    for k in xrange(1, robo.NL):
        param_vec = robo.get_inert_param(k)
        F = ParamsInit.init_vec(robo)
        N = ParamsInit.init_vec(robo)
        for i in xrange(10):
            if param_vec[i] == tools.ZERO:
                continue
            # change link names according to current non-zero parameter
            robo_tmp.num = [str(l) + str(param_vec[i])
                            for l in xrange(k + 1)]
            # set the parameter to 1
            mask = sympy.zeros(10, 1)
            mask[i] = 1
            robo_tmp.put_inert_param(mask, k)
            # compute the total forcec of the link k
            compute_wrench(robo_tmp, symo, k, w, wdot, U, vdot, F, N)
            # init external forces
            Fex = copy(robo.Fex)
            Nex = copy(robo.Nex)
            for j in reversed(xrange(k + 1)):
                compute_joint_wrench(robo_tmp, symo, j, antRj, antPj,
                                     vdot, Fjnt, Njnt, F, N, Fex, Nex)
            for j in xrange(k + 1):
                compute_torque(robo_tmp, symo, j, Fjnt, Njnt, 'DG')
        # reset all the parameters to zero
        robo_tmp.put_inert_param(sympy.zeros(10, 1), k)
        # compute model for the joint parameters
        compute_joint_torque_deriv(symo, robo.IA[k],
                                   robo.qddot[k], k)
        compute_joint_torque_deriv(symo, robo.FS[k],
                                   sympy.sign(robo.qdot[k]), k)
        compute_joint_torque_deriv(symo, robo.FV[k],
                                   robo.qdot[k], k)
    # closing the output file
    symo.file_close()
    return symo


def direct_dynamic_NE(robo):
    """Computes Direct Dynamic Model using
    Newton-Euler formulation

    Parameters
    ==========
    robo : Robot
        Instance of robot description container

    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    wi = ParamsInit.init_vec(robo)
        # antecedent angular velocity, projected into jth frame
    w = ParamsInit.init_w(robo)
    jaj = ParamsInit.init_vec(robo, 6)
    jTant = ParamsInit.init_mat(robo, 6)   # Twist transform list of Matrices 6x6
    beta_star = ParamsInit.init_vec(robo, 6)
    grandJ = ParamsInit.init_mat(robo, 6)
    link_acc = ParamsInit.init_vec(robo, 6)
    H_inv = ParamsInit.init_scalar(robo)
    juj = ParamsInit.init_vec(robo, 6)   # Jj*aj / Hj
    Tau = ParamsInit.init_scalar(robo)
    grandVp = ParamsInit.init_vec(robo, 6)
    grandVp.append(Matrix([robo.vdot0 - robo.G, robo.w0]))
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'ddm')
    title = 'Direct dynamic model using Newton - Euler Algorith'
    symo.write_params_table(robo, title, inert=True, dynam=True)

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
        if robo.ant[j] != - 1:
            compute_beta_J_star(robo, symo, j, grandJ, jaj, juj, Tau,
                                beta_star, jTant, link_acc)
    for j in xrange(1, robo.NL):
        compute_acceleration(robo, symo, j, jTant, grandVp,
                             juj, H_inv, jaj, Tau, link_acc)
    for j in xrange(1, robo.NL):
        compute_coupled_forces(robo, symo, j, grandVp, grandJ, beta_star)
    symo.file_close()
    return symo


def inertia_matrix(robo):
    """Computes Inertia Matrix using composed link

    Parameters
    ==========
    robo : Robot
        Instance of robot description container

    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    Jplus, MSplus, Mplus = ParamsInit.init_jplus(robo)
    AJE1 = ParamsInit.init_vec(robo)
    f = ParamsInit.init_vec(robo, ext=1)
    n = ParamsInit.init_vec(robo, ext=1)
    A = sympy.zeros(robo.NL, robo.NL)
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'inm')
    title = 'Inertia Matrix using composite links'
    symo.write_params_table(robo, title, inert=True, dynam=True)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    for j in reversed(xrange(-1, robo.NL)):
        replace_Jplus(robo, symo, j, Jplus, MSplus, Mplus)
        if j != - 1:
            compute_Jplus(robo, symo, j, antRj, antPj,
                          Jplus, MSplus, Mplus, AJE1)
    for j in xrange(1, robo.NL):
        compute_A_diagonal(robo, symo, j, Jplus, MSplus, Mplus, f, n, A)
        ka = j
        while ka != - 1:
            k = ka
            ka = robo.ant[ka]
            compute_A_triangle(robo, symo, j, k, ka,
                               antRj, antPj, f, n, A, AJE1)
    symo.mat_replace(A, 'A', forced=True, symmet=True)
    J_base = inertia_spatial(Jplus[-1], MSplus[-1], Mplus[-1])
    symo.mat_replace(J_base, 'JP', 0, forced=True, symmet=True)
    symo.file_close()
    return symo


def inverse_dynamic_NE(robo):
    """Computes Inverse Dynamic Model using
    Newton-Euler formulation

    Parameters
    ==========
    robo : Robot
        Instance of robot description container

    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'idm')
    title = 'Inverse dynamic model using Newton - Euler Algorith'
    symo.write_params_table(robo, title, inert=True, dynam=True)
    Newton_Euler(robo, symo)
    symo.file_close()
    return symo


def pseudo_force_NE(robo):
    """Computes Coriolis, Centrifugal, Gravity, Friction and external
    torques using Newton-Euler formulation

    Parameters
    ==========
    robo : Robot
        Instance of robot description container

    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    robo_pseudo = deepcopy(robo)
    robo_pseudo.qddot = sympy.zeros(robo_pseudo.NL, 1)
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'ccg')
    title = 'Pseudo forces using Newton - Euler Algorith'
    symo.write_params_table(robo, title, inert=True, dynam=True)
    Newton_Euler(robo_pseudo, symo)
    symo.file_close()
    return symo


def compute_wrench(robo, symo, j, w, wdot, U, vdot, F, N):
    """Internal function. Computes total wrench (torques and forces)
    of the link j

    Notes
    =====
    F, N are the output parameters
    """
    F[j] = robo.M[j]*vdot[j] + U[j]*robo.MS[j]
    symo.mat_replace(F[j], 'F', j)
    Psi = symo.mat_replace(robo.J[j]*w[j], 'PSI', j)
    N[j] = robo.J[j]*wdot[j] + tools.skew(w[j])*Psi
    symo.mat_replace(N[j], 'No', j)


def compute_joint_wrench(robo, symo, j, antRj, antPj, vdot,
                         Fjnt, Njnt, F, N, Fex, Nex):
    """Internal function. Computes wrench (torques and forces)
    of the joint j

    Notes
    =====
    Fjnt, Njnt, Fex, Nex are the output parameters
    """
    Fjnt[j] = symo.mat_replace(F[j] + Fex[j], 'E', j)
    Njnt[j] = N[j] + Nex[j] + tools.skew(robo.MS[j])*vdot[j]
    symo.mat_replace(Njnt[j], 'N', j)
    f_ant = symo.mat_replace(antRj[j]*Fjnt[j], 'FDI', j)
    if robo.ant[j] != - 1:
        Fex[robo.ant[j]] += f_ant
        Nex[robo.ant[j]] += antRj[j]*Njnt[j] + tools.skew(antPj[j])*f_ant


def compute_torque(robo, symo, j, Fjnt, Njnt, name='GAM'):
    """Internal function. Computes actuation torques - projection of
    joint wrench on the joint axis
    """
    if robo.sigma[j] != 2:
        tau = (robo.sigma[j]*Fjnt[j] + (1 - robo.sigma[j])*Njnt[j])
        tau_total = tau[2] + robo.fric_s(j) + robo.fric_v(j) + robo.tau_ia(j)
        symo.replace(tau_total, name, j, forced=True)


def inertia_spatial(J, MS, M):
    return Matrix([(M*sympy.eye(3)).row_join(tools.skew(MS).T), tools.skew(MS).row_join(J)])


def compute_joint_torque_deriv(symo, param, arg, index):
    """Internal function. Computes joint reactive torques
    in case if the parameter is 1

    Parameters
    ==========
    symo : symbolmgr.SymbolManager
        symbol manager
    param : var
        Dynamic parameter
    arg : var
        The real torque is equal to arg*param
    index : strig
        identifies the parameter in the sybstituted symbol's name
    """
    if param != tools.ZERO and arg != tools.ZERO:
        index = str(index) + str(param)
        symo.replace(arg, 'DG', index, forced=True)


def compute_beta(robo, symo, j, w, beta_star):
    """Internal function. Computes link's wrench when
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
    """Internal function. Computes link's accelerations when
    the joint accelerations are zero

    Notes
    =====
    link_acc is the output parameter
    """
    E1 = symo.mat_replace(tools.skew(wi[j])*Matrix([0, 0, robo.qdot[j]]),
                          'WQ', j)
    E2 = (1 - robo.sigma[j])*E1
    E3 = 2*robo.sigma[j]*E1
    E4 = tools.skew(w[robo.ant[j]])*antPj[j]
    E5 = tools.skew(w[robo.ant[j]])*E4
    E6 = antRj[j].T*E5
    E7 = symo.mat_replace(E6 + E3, 'LW', j)
    link_acc[j] = Matrix([E7, E2])


def replace_beta_J_star(robo, symo, j, grandJ, beta_star):
    """Internal function. Makes symbol substitution in beta_star
    and grandJ
    """
    grandJ[j] = symo.mat_replace(grandJ[j], 'MJE', j, symmet=True)
    beta_star[j] = symo.mat_replace(beta_star[j], 'VBE', j)


def compute_Tau(robo, symo, j, grandJ, beta_star, jaj, juj, H_inv, Tau):
    """Internal function. Computes intermediat dynamic variables

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
    """Internal function. Computes intermediat dynamic variables

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
    """Internal function. Computes joint accelerations and links' twists

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
    qddot = symo.replace(qddot, "QDP", j, forced=True)
    grandVp[j] = (grandR + qddot*jaj[j])
    grandVp[j][3:, 0] = symo.mat_replace(grandVp[j][3:, 0], 'WP', j)
    grandVp[j][:3, 0] = symo.mat_replace(grandVp[j][:3, 0], 'VP', j)


def compute_coupled_forces(robo, symo, j, grandVp, grandJ, beta_star):
    """Internal function.
    """
    E3 = symo.mat_replace(grandJ[j]*grandVp[j], 'DY', j)
    couplforce = E3 - beta_star[j]
    symo.mat_replace(couplforce[3:, 0], 'N', j)
    symo.mat_replace(couplforce[:3, 0], 'E', j)


def replace_Jplus(robo, symo, j, Jplus, MSplus, Mplus):
    """Internal function. Makes symbol substitutions inertia parameters
    """
    symo.mat_replace(Jplus[j], 'JP', j)
    symo.mat_replace(MSplus[j], 'MSP', j)
    Mplus[j] = symo.replace(Mplus[j], 'MP', j)


def compute_Jplus(robo, symo, j, antRj, antPj, Jplus, MSplus, Mplus, AJE1):
    """Internal function. Computes inertia parameters of composed link

    Notes
    =====
    Jplus, MSplus, Mplus are the output parameters
    """
    hat_antPj = tools.skew(antPj[j])
    antMSj = symo.mat_replace(antRj[j]*MSplus[j], 'AS', j)
    E1 = symo.mat_replace(antRj[j]*Jplus[j], 'AJ', j)
    AJE1[j] = E1[:, 2]
    E2 = symo.mat_replace(E1*antRj[j].T, 'AJA', j)
    E3 = symo.mat_replace(hat_antPj*tools.skew(antMSj), 'PAS', j)
    Jplus[robo.ant[j]] += E2 - (E3 + E3.T) + hat_antPj*hat_antPj.T*Mplus[j]
    MSplus[robo.ant[j]] += antMSj + antPj[j]*Mplus[j]
    Mplus[robo.ant[j]] += Mplus[j]


def compute_A_diagonal(robo, symo, j, Jplus, MSplus, Mplus, f, n, A):
    """Internal function. Computes diagonal elements
    of the inertia matrix

    Notes
    =====
    f, n, A are the output parameters
    """
    if robo.sigma[j] == 0:
        f[j] = Matrix([-MSplus[j][1], MSplus[j][0], 0])
        n[j] = Jplus[j][:, 2]
        A[j, j] = Jplus[j][2, 2] + robo.IA[j]
    elif robo.sigma[j] == 1:
        f[j] = Matrix([0, 0, Mplus[j]])
        n[j] = Matrix([MSplus[j][1], - MSplus[j][0], 0])
        A[j, j] = Mplus[j] + robo.IA[j]
    symo.mat_replace(f[j], 'E' + chars[j], j)
    symo.mat_replace(n[j], 'N' + chars[j], j)


def compute_A_triangle(robo, symo, j, k, ka, antRj, antPj, f, n, A, AJE1):
    """Internal function. Computes elements below and above diagonal
    of the inertia matrix

    Notes
    =====
    f, n, A are the output parameters
    """
    f[ka] = antRj[k]*f[k]
    if k == j and robo.sigma[j] == 0:
        n[ka] = AJE1[k] + tools.skew(antPj[k])*f[k]
    else:
        n[ka] = antRj[k]*n[k] + tools.skew(antPj[k])*f[k]
    if ka == - 1:
        symo.mat_replace(f[ka], 'AV0')
        symo.mat_replace(n[ka], 'AW0')
    else:
        symo.mat_replace(f[ka], 'E' + chars[j], ka)
        symo.mat_replace(n[ka], 'N' + chars[j], ka)
        if robo.sigma[ka] == 0:
            A[j, ka] = n[ka][2]
        elif robo.sigma[ka] == 1:
            A[j, ka] = f[ka][2]
        A[ka, j] = A[j, ka]


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
            compute_lambda(robo, symo, j, antRj, antPj)
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


