#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module contains the labels used in the user interface as a dict.
The main purpose of this to represent the labels symbolically and use it
in multiple places. Also this gives the advantage of modifying and
maintaining the labels easily since all the labels are placed in the 
same place.
"""


MAIN_WIN = dict(
    prog_name = "OpenSYMORO",
    window_title = "OpenSYMORO - SYmbolic MOdelling of RObots"
)

MAIN_MENU = dict(
    file_menu = "&File",
    geom_menu = "&Geometric",
    kin_menu = "&Kinematic",
    dyn_menu = "&Dynamic",
    iden_menu = "&Identification",
    optim_menu = "&Optimiser",
    viz_menu = "&Visualisation"
)

VIZ_MENU = dict(m_viz = "Visualisation")

IDEN_MENU = dict(
    m_base_inertial_params = "Base Inertial parameters",
    m_dyn_iden_model = "Dynamic Identification Model",
    m_energy_iden_model = "Energy Identification Model"
)

DYN_MENU = dict(
    m_idym = "Inverse Dynamic Model",
    m_inertia_matrix = "Inertia matrix",
    m_h_term = "Centrifugal, Coriolis & Gravity torques",
    m_ddym = "Direct Dynamic Model"
)

KIN_MENU = dict(
    m_jac_matrix = "Jacobian matrix",
    m_determinant = "Determinant of a Jacobian",
    m_kin_constraint = "Kinematic constraint equation of loops",
    m_vel = "Velocities",
    m_acc = "Accelerations",
    m_jpqp = "Jpqp"
)

GEOM_MENU = dict(
    m_trans_matrix = "Transformation matrix",
    m_fast_dgm = "Fast Geometric model",
    m_igm_paul = "IGM - Paul method",
    m_geom_constraint = "Geometric constraint equation of loops"
)

FILE_MENU = dict(
    m_new = "&New",
    m_open = "&Open",
    m_save = "&Save",
    m_save_as = "Save &As",
    m_pref = "Preferences\t(unavailable)",
    m_exit = "E&xit"
)

