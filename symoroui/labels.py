#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module contains the labels used in the user interface as a dict.
The main purpose of this to represent the labels symbolically and use it
in multiple places. Also this gives the advantage of modifying and
maintaining the labels easily since all the labels are placed in the 
same place.
"""


from collections import namedtuple
from collections import OrderedDict


# main window
MAIN_WIN = dict(
    prog_name = "SYMORO",
    window_title = "SYMORO - SYmbolic MOdelling of RObots"
)

# interface contents
# TODO: add event handlers to dict entries as well
BOX_TITLES = dict(
    robot_des = "Robot Description",
    robot_type = "Robot Type",
    gravity = "Gravity Components",
    location = "Robot Location",
    geom_params = "Geometric Parameters",
    dyn_params = "Dynamic Parameters and External Forces",
    base_vel_acc = "Velocity and Acceleration of the base",
    joint_vel_acc = "Joint Velocity and Acceleration"
)
# named tuple to hold the content field entries
FieldEntry = namedtuple('FieldEntry', ['label', 'name', 'control'])
# joint velocity and acceleration params
JOINT_VEL_ACC = OrderedDict([
    ('joint', FieldEntry('Joint', 'joint', 'cmb')),
    ('qp', FieldEntry('QP', 'QP', 'txt')),
    ('qdp', FieldEntry('QDP', 'QDP', 'txt')),
    ('gam', FieldEntry('GAM', 'GAM', 'txt')),
])
# base velocity and acceleration params
BASE_VEL_ACC = OrderedDict([
    ('wx', FieldEntry('W0X', 'W0X', 'txt')),
    ('wy', FieldEntry('W0Y', 'W0Y', 'txt')),
    ('wz', FieldEntry('W0Z', 'W0Z', 'txt')),
    ('wpx', FieldEntry('WP0X', 'WP0X', 'txt')),
    ('wpy', FieldEntry('WP0Y', 'WP0Y', 'txt')),
    ('wpz', FieldEntry('WP0Z', 'WP0Z', 'txt')),
    ('vx', FieldEntry('V0X', 'V0X', 'txt')),
    ('vy', FieldEntry('V0Y', 'V0Y', 'txt')),
    ('vz', FieldEntry('V0Z', 'V0Z', 'txt')),
    ('vpx', FieldEntry('VP0X', 'VP0X', 'txt')),
    ('vpy', FieldEntry('VP0Y', 'VP0Y', 'txt')),
    ('vpz', FieldEntry('VP0Z', 'VP0Z', 'txt'))
])
# inertial params
DYN_PARAMS_I = OrderedDict([
    ('xx', FieldEntry('XX', 'XX', 'txt')),
    ('xy', FieldEntry('XY', 'XY', 'txt')),
    ('xz', FieldEntry('XZ', 'XZ', 'txt')),
    ('yy', FieldEntry('YY', 'YY', 'txt')),
    ('yz', FieldEntry('YZ', 'YZ', 'txt')),
    ('zz', FieldEntry('ZZ', 'ZZ', 'txt'))
])
# mass tensor params
DYN_PARAMS_M = OrderedDict([
    ('mx', FieldEntry('MX', 'MX', 'txt')),
    ('my', FieldEntry('MY', 'MY', 'txt')),
    ('mz', FieldEntry('MZ', 'MZ', 'txt')),
    ('m', FieldEntry('M', 'M', 'txt'))
])
# friction and rotor inertia params
DYN_PARAMS_X = OrderedDict([
    ('ia', FieldEntry('IA', 'IA', 'txt')),
    ('fc', FieldEntry('FS', 'FS', 'txt')),
    ('fv', FieldEntry('FV', 'FV', 'txt'))
])
# external force, moments params
DYN_PARAMS_F = OrderedDict([
    ('ex_fx', FieldEntry('FX', 'FX', 'txt')),
    ('ex_fy', FieldEntry('FY', 'FY', 'txt')),
    ('ex_fz', FieldEntry('FZ', 'FZ', 'txt')),
    ('ex_mx', FieldEntry('CX', 'CX', 'txt')),
    ('ex_my', FieldEntry('CY', 'CY', 'txt')),
    ('ex_mz', FieldEntry('CZ', 'CZ', 'txt'))
])
# dynamic params got by concatenation
DYN_PARAMS = OrderedDict(
    [('link', FieldEntry('Link', 'link', 'cmb'))] + \
    DYN_PARAMS_I.items() + \
    DYN_PARAMS_M.items() + \
    DYN_PARAMS_X.items() + \
    DYN_PARAMS_F.items()
)
# geometric params
GEOM_PARAMS = OrderedDict([
    ('frame', FieldEntry('Frame', 'frame', 'cmb')),
    ('ant', FieldEntry('ant', 'ant', 'cmb')),
    ('sigma', FieldEntry('sigma', 'sigma', 'cmb')),
    ('mu', FieldEntry('mu', 'mu', 'cmb')),
    ('gamma', FieldEntry('gamma', 'gamma', 'txt')),
    ('b', FieldEntry('b', 'b', 'txt')),
    ('alpha', FieldEntry('alpha', 'alpha', 'txt')),
    ('d', FieldEntry('d', 'd', 'txt')),
    ('theta', FieldEntry('theta', 'theta', 'txt')),
    ('r', FieldEntry('r', 'r', 'txt'))
])
# gravity component params
GRAVITY_CMPNTS = OrderedDict([
    ('gx', FieldEntry('GX', 'GX', 'txt')),
    ('gy', FieldEntry('GY', 'GY', 'txt')),
    ('gz', FieldEntry('GZ', 'GZ', 'txt'))
])
# robot type params
ROBOT_TYPE = OrderedDict([
    ('name', FieldEntry('Name of the robot:', 'name', 'lbl')),
    ('num_links', FieldEntry('Number of moving links:', 'NL', 'lbl')),
    ('num_joints', FieldEntry('Number of joints:', 'NJ', 'lbl')), 
    ('num_frames', FieldEntry('Number of frames:', 'NF', 'lbl')),
    ('structure', FieldEntry('Type of structure:', 'type', 'lbl')),
    ('is_mobile', FieldEntry('Is Mobile:', 'mobile', 'lbl')),
    ('num_loops', FieldEntry('Number of closed loops:', 'loops', 'lbl'))
])

# menu bar
# TODO: add event handlers to dict entries as well
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
    m_pref = "Preferences -- (unavailable)",
    m_exit = "E&xit"
)


