# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module contains the functions for the computation of Dynamic
Identification model.
"""


import copy

import sympy

from pysymoro.geometry import compute_rot_trans
from pysymoro.kinematics import compute_vel_acc
from symoroutils.paramsinit import ParamsInit
from symoroutils import tools


def get_symbol(symbol, name, index, element=''):
    return symbol + name.format(index=index, element=element)


def vec_replace_wrapper(symo, vec, symbol, name, index, forced=False):
    new_vec = sympy.zeros(vec.rows, 1)
    for idx in xrange(vec.rows):
        sym = get_symbol(symbol, name, index, idx+1)
        new_vec[idx] = symo.replace(vec[idx], sym, forced=forced)
    return new_vec


def _compute_dynamic_wrench(robo, symo, name, j, w, wdot, U, vdot, F, N):
    """
    Compute total wrench of link j (internal function).

    Note:
        F, N are the output parameters
    """
    F[j] = (robo.M[j] * vdot[j]) + (U[j] * robo.MS[j])
    F[j] = vec_replace_wrapper(symo, F[j], 'F', name, j)
    Psi = robo.J[j] * w[j]
    Psi = vec_replace_wrapper(symo, Psi, 'PSI', name, j)
    N[j] = (robo.J[j] * wdot[j]) + (tools.skew(w[j]) * Psi)
    N[j] = vec_replace_wrapper(symo, N[j], 'No', name, j)


def _compute_reaction_wrench(
    robo, symo, name, j, antRj, antPj, vdot, F, N, Fjnt, Njnt, Fex, Nex
):
    """
    Compute reaction wrench (for default Newton-Euler) of joint j
    (internal function).

    Note:
        Fjnt, Njnt, Fex, Nex are the output parameters
    """
    i = robo.ant[j]
    Fjnt[j] = F[j] + Fex[j]
    Fjnt[j] = vec_replace_wrapper(symo, Fjnt[j], 'E', name, j)
    Njnt[j] = N[j] + Nex[j] + (tools.skew(robo.MS[j]) * vdot[j])
    Njnt[j] = vec_replace_wrapper(symo, Njnt[j], 'N', name, j)
    f_ant = antRj[j] * Fjnt[j]
    f_ant = vec_replace_wrapper(symo, f_ant, 'FDI', name, j)
    if i != -1:
        Fex[i] = Fex[i] + f_ant
        Nex[i] = Nex[i] + \
            (antRj[j] * Njnt[j]) + (tools.skew(antPj[j]) * f_ant)


def _compute_base_reaction_wrench(
    robo, symo, name, antRj, antPj, vdot, F, N, Fex, Nex, Fjnt, Njnt
):
    """
    Compute reaction wrench (for default Newton-Euler) on the base
    (internal function).

    Note:
        Fjnt, Njnt are the output parameters
    """
    j = 0
    Fjnt[j] = F[j] + Fex[j]
    Fjnt[j] = vec_replace_wrapper(
        symo, Fjnt[j], 'DE', name, j, forced=True
    )
    Njnt[j] = N[j] + Nex[j] + (tools.skew(robo.MS[j]) * vdot[j])
    Njnt[j] = vec_replace_wrapper(
        symo, Njnt[j], 'DN', name, j, forced=True
    )


def _compute_joint_torque(robo, symo, name, j, Fjnt, Njnt):
    """
    Compute actuator torques - projection of joint wrench on the joint
    axis (internal function).
    """
    if robo.sigma[j] == 2:
        tau_total = 0
    else:
        tau = (robo.sigma[j] * Fjnt[j]) + ((1 - robo.sigma[j]) * Njnt[j])
        fric_rotor = robo.fric_s(j) + robo.fric_v(j) + robo.tau_ia(j)
        tau_total = tau[2] + fric_rotor
    symo.replace(tau_total, get_symbol('DG', name, j), forced=True)


def _compute_joint_torque_deriv(symo, param, arg, index):
    """Compute joint reactive torque if the parameter is 1

    Parameters:
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


def dynamic_identification_model(robo, symo):
    """
    Compute the Dynamic Identification model of a robot using
    Newton-Euler algorithm.
    """
    # init forces vectors
    Fjnt = ParamsInit.init_vec(robo)
    Njnt = ParamsInit.init_vec(robo)
    # init transformation
    antRj, antPj = compute_rot_trans(robo, symo)
    # init velocities and accelerations
    w, wdot, vdot, U = compute_vel_acc(
        robo, symo, antRj, antPj, floating=True, gravity=True
    )
    # virtual robot with only one non-zero parameter at once
    robo_tmp = copy.deepcopy(robo)
    robo_tmp.IA = sympy.zeros(robo.NL, 1)
    robo_tmp.FV = sympy.zeros(robo.NL, 1)
    robo_tmp.FS = sympy.zeros(robo.NL, 1)
    # start link number
    is_fixed = False if robo.is_floating or robo.is_mobile else True
    start_link = 0
    for k in xrange(start_link, robo.NL):
        param_vec = robo.get_inert_param(k)
        F = ParamsInit.init_vec(robo)
        N = ParamsInit.init_vec(robo)
        for i in xrange(10):
            if param_vec[i] == tools.ZERO:
                continue
            # change link names according to current non-zero parameter
            name = '{index}{element}' + str(param_vec[i])
            # set the parameter to 1
            mask = sympy.zeros(10, 1)
            mask[i] = 1
            robo_tmp.put_inert_param(mask, k)
            # compute the total forcec of the link k
            _compute_dynamic_wrench(
                robo_tmp, symo, name, k, w, wdot, U, vdot, F, N
            )
            # init external forces
            Fex = ParamsInit.init_vec(robo)
            Nex = ParamsInit.init_vec(robo)
            for j in reversed(xrange(1, k + 1)):
                _compute_reaction_wrench(
                    robo_tmp, symo, name, j, antRj, antPj,
                    vdot, F, N, Fjnt, Njnt, Fex, Nex
                )
            # reaction wrench for base
            _compute_base_reaction_wrench(
                robo_tmp, symo, name, antRj,antPj,
                vdot, F, N, Fex, Nex, Fjnt, Njnt
            )
            for j in xrange(1, k + 1):
                _compute_joint_torque(robo_tmp, symo, name, j, Fjnt, Njnt)
        # reset all the parameters to zero
        robo_tmp.put_inert_param(sympy.zeros(10, 1), k)
        # compute model for the joint parameters
        # avoid these parameters for link 0
        if k == 0: continue
        _compute_joint_torque_deriv(
            symo, robo.IA[k], robo.qddot[k], k
        )
        _compute_joint_torque_deriv(
            symo, robo.FS[k], sympy.sign(robo.qdot[k]), k
        )
        _compute_joint_torque_deriv(
            symo, robo.FV[k], robo.qdot[k], k
        )
    return symo


