# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module creates the dialog boxes to specify parameters for
dynamic model calculations.
"""
import copy

from sympy import zeros

from pysymoro import nealgos, inertia, baseparams, dyniden
from symoroutils import filemgr, symbolmgr


def compute_idym(robo):
    """
    Compute the Inverse Dynamic Model of the robot using the
    recursive Newton-Euler algorithm. Also choose the Newton-Euler
    algorithm based on the robot type.
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'idm')
    title = "Inverse Dynamic Model using Newton-Euler Algorithm\n"
    if 1 in robo.eta:
        # with flexible joints
        title = title + "Robot with flexible joints\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.flexible_inverse_dynmodel(robo, symo)
    elif robo.is_floating:
        # with rigid joints and floating base
        title = title + "Robot with rigid joints and floating base\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.composite_inverse_dynmodel(robo, symo)
    elif robo.is_mobile:
        # mobile robot with rigid joints - known base acceleration
        title = title + "Robot with mobile base (Vdot0 is known)\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.mobile_inverse_dynmodel(robo, symo)
    else:
        # with rigid joints and fixed base
        title = title + "Robot with rigid joints and fixed base\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.fixed_inverse_dynmodel(robo, symo)
    symo.file_close()
    return symo


def compute_inertiamatrix(robo):
    """
    Compute the Inertia Matrix of the robot using the Composite link
    algorithm.
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'inm')
    title = "Inertia matrix using Composite links algorithm\n"
    if robo.is_floating or robo.is_mobile:
        # with floating or mobile base
        title = title + "Robot with floating/mobile base\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        inertia.floating_inertia_matrix(robo, symo)
    else:
        # with fixed base
        title = title + "Robot with fixed base\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        inertia.fixed_inertia_matrix(robo, symo)
    symo.file_close()
    return symo


def compute_ddym(robo):
    """
    Compute the Direct Dynamic Model of the robot using the
    recursive Newton-Euler algorithm.
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'ddm')
    title = "Direct Dynamic Model using Newton-Euler Algorithm\n"
    if robo.is_floating:
        # with floating base
        title = title + "Robot with floating base\n"
    else:
        # with fixed base
        title = title + "Robot with fixed base\n"
    symo.write_params_table(robo, title, inert=True, dynam=True)
    nealgos.direct_dynmodel(robo, symo)
    symo.file_close()
    return symo


def compute_pseudotorques(robo):
    """
    Compute Coriolis, Centrifugal, Gravity, Friction and external
    torques using Newton-Euler algortihm.
    """
    pseudo_robo = copy.deepcopy(robo)
    pseudo_robo.qddot = zeros(pseudo_robo.NL, 1)
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'ccg')
    title = "Pseudo forces using Newton-Euler Algorithm\n"
    if 1 in pseudo_robo.eta:
        # with flexible joints
        title = title + "Robot with flexible joints\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.flexible_inverse_dynmodel(pseudo_robo, symo)
    elif pseudo_robo.is_floating:
        # with rigid joints and floating base
        title = title + "Robot with rigid joints and floating base\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.composite_inverse_dynmodel(pseudo_robo, symo)
    elif pseudo_robo.is_mobile:
        # mobile robot with rigid joints - known base acceleration
        title = title + "Robot with mobile base (Vdot0 is known)\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.mobile_inverse_dynmodel(pseudo_robo, symo)
    else:
        # with rigid joints and fixed base
        title = title + "Robot with rigid joints and fixed base\n"
        symo.write_params_table(robo, title, inert=True, dynam=True)
        nealgos.fixed_inverse_dynmodel(pseudo_robo, symo)
    symo.file_close()
    return symo


def compute_baseparams(robo):
    """
    Compute the Base Inertial Parameters of the robot.
    """
    base_robo = copy.deepcopy(robo)
    symo = symbolmgr.SymbolManager()
    symo.file_open(base_robo, 'regp')
    title = "Base Inertial Parameters equations"
    symo.write_params_table(
        base_robo, title, inert=True, dynam=True
    )
    # compute base inertial params
    baseparams.base_inertial_parameters(base_robo, symo)
    symo.write_line()
    title = "Grouped inertial parameters"
    symo.write_params_table(
        base_robo, title, inert=True, equations=False
    )
    symo.file_close()
    # set new name for robot with base params
    base_robo.name = base_robo.name + "_base"
    file_path = filemgr.get_file_path(base_robo)
    base_robo.set_par_file_path(file_path)
    return symo, base_robo


def compute_dynidenmodel(robo):
    """
    Compute the Dynamic Identification model of the robot.
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'dim')
    title = "Dynamic Identification Model (Newton-Euler method)"
    symo.write_params_table(robo, title, inert=True, dynam=True)
    dyniden.dynamic_identification_model(robo, symo)
    symo.file_close()
    return symo
