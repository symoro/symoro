# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


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
    prog_name='SYMORO',
    window_title='SYMORO - SYmbolic MOdelling of RObots'
)

# interface contents
# TODO: add event handlers to dict entries as well
BOX_TITLES = dict(
    robot_des='Robot Description',
    robot_type='Robot Type',
    gravity='Gravity Components',
    location='Robot Location',
    geom_params='Geometric Parameters',
    dyn_params='Dynamic Parameters and External Forces',
    base_vel_acc='Velocity and Acceleration of the base',
    joint_params='Joint Stiffness, Velocity and Acceleration'
)
# named tuple to hold the content field entries
FieldEntry = namedtuple(
    'FieldEntry', ['label', 'name', 'control', 'place', 'handler', 'id']
)
# joint velocity and acceleration params
JOINT_PARAMS = OrderedDict([
    ('joint', FieldEntry('Joint', 'joint', 'cmb', (0, 0), 'OnJointChanged', -1)),
    ('eta', FieldEntry('eta', 'eta', 'cmb', (0, 1), 'OnSpeedChanged', -1)),
    ('stiff', FieldEntry('k', 'k', 'txt', (1, 0), 'OnSpeedChanged', -1)),
    ('gam', FieldEntry('GAM', 'GAM', 'txt', (1, 1), 'OnSpeedChanged', -1)),
    ('qp', FieldEntry('QP', 'QP', 'txt', (2, 0), 'OnSpeedChanged', -1)),
    ('qdp', FieldEntry('QDP', 'QDP', 'txt', (2, 1), 'OnSpeedChanged', -1))
])
# base velocity and acceleration params
BASE_VEL_ACC = OrderedDict([
    ('vx', FieldEntry('VXb', 'V0X', 'txt', (0, 0), 'OnBaseTwistChanged', 0)),
    ('vy', FieldEntry('VYb', 'V0Y', 'txt', (1, 0), 'OnBaseTwistChanged', 1)),
    ('vz', FieldEntry('VZb', 'V0Z', 'txt', (2, 0), 'OnBaseTwistChanged', 2)),
    ('wx', FieldEntry('WXb', 'W0X', 'txt', (0, 1), 'OnBaseTwistChanged', 0)),
    ('wy', FieldEntry('WYb', 'W0Y', 'txt', (1, 1), 'OnBaseTwistChanged', 1)),
    ('wz', FieldEntry('WZb', 'W0Z', 'txt', (2, 1), 'OnBaseTwistChanged', 2)),
    ('vpx', FieldEntry('VPXb', 'VP0X', 'txt', (0, 2), 'OnBaseTwistChanged', 0)),
    ('vpy', FieldEntry('VPYb', 'VP0Y', 'txt', (1, 2), 'OnBaseTwistChanged', 1)),
    ('vpz', FieldEntry('VPZb', 'VP0Z', 'txt', (2, 2), 'OnBaseTwistChanged', 2)),
    ('wpx', FieldEntry('WPXb', 'WP0X', 'txt', (0, 3), 'OnBaseTwistChanged', 0)),
    ('wpy', FieldEntry('WPYb', 'WP0Y', 'txt', (1, 3), 'OnBaseTwistChanged', 1)),
    ('wpz', FieldEntry('WPZb', 'WP0Z', 'txt', (2, 3), 'OnBaseTwistChanged', 2))
])
# inertial params
DYN_PARAMS_I = OrderedDict([
    ('xx', FieldEntry('XX', 'XX', 'txt', (0, 0), 'OnDynParamChanged', -1)),
    ('xy', FieldEntry('XY', 'XY', 'txt', (0, 1), 'OnDynParamChanged', -1)),
    ('xz', FieldEntry('XZ', 'XZ', 'txt', (0, 2), 'OnDynParamChanged', -1)),
    ('yy', FieldEntry('YY', 'YY', 'txt', (0, 3), 'OnDynParamChanged', -1)),
    ('yz', FieldEntry('YZ', 'YZ', 'txt', (0, 4), 'OnDynParamChanged', -1)),
    ('zz', FieldEntry('ZZ', 'ZZ', 'txt', (0, 5), 'OnDynParamChanged', -1))
])
# mass tensor params
DYN_PARAMS_M = OrderedDict([
    ('msx', FieldEntry('MX', 'MX', 'txt', (1, 0), 'OnDynParamChanged', -1)),
    ('msy', FieldEntry('MY', 'MY', 'txt', (1, 1), 'OnDynParamChanged', -1)),
    ('msz', FieldEntry('MZ', 'MZ', 'txt', (1, 2), 'OnDynParamChanged', -1)),
    ('mass', FieldEntry('M', 'M', 'txt', (1, 3), 'OnDynParamChanged', -1))
])
# friction and rotor inertia params
DYN_PARAMS_X = OrderedDict([
    ('ia', FieldEntry('IA', 'IA', 'txt', (1, 4), 'OnDynParamChanged', -1)),
    ('frc', FieldEntry('FS', 'FS', 'txt', (1, 5), 'OnDynParamChanged', -1)),
    ('frv', FieldEntry('FV', 'FV', 'txt', (1, 6), 'OnDynParamChanged', -1))
])
# external force, moments params
DYN_PARAMS_F = OrderedDict([
    ('fx_ext', FieldEntry('FX', 'FX', 'txt', (2, 0), 'OnDynParamChanged', -1)),
    ('fy_ext', FieldEntry('FY', 'FY', 'txt', (2, 1), 'OnDynParamChanged', -1)),
    ('fz_ext', FieldEntry('FZ', 'FZ', 'txt', (2, 2), 'OnDynParamChanged', -1)),
    ('mx_ext', FieldEntry('CX', 'CX', 'txt', (2, 3), 'OnDynParamChanged', -1)),
    ('my_ext', FieldEntry('CY', 'CY', 'txt', (2, 4), 'OnDynParamChanged', -1)),
    ('mz_ext', FieldEntry('CZ', 'CZ', 'txt', (2, 5), 'OnDynParamChanged', -1))
])
# dynamic params got by concatenation
DYN_PARAMS = OrderedDict(
    [('link', FieldEntry('Link', 'link', 'cmb', None, 'OnLinkChanged', -1))])
DYN_PARAMS.update(DYN_PARAMS_I)
DYN_PARAMS.update(DYN_PARAMS_M)
DYN_PARAMS.update(DYN_PARAMS_X)
DYN_PARAMS.update(DYN_PARAMS_F)
# geometric params
GEOM_PARAMS = OrderedDict([
    ('frame', FieldEntry('Frame', 'frame', 'cmb', (0, 0), 'OnFrameChanged', -1)),
    ('ant', FieldEntry('ant', 'ant', 'cmb', (1, 0), 'OnGeoParamChanged', -1)),
    ('sigma', FieldEntry('sigma', 'sigma', 'cmb', (0, 1), 'OnGeoParamChanged', -1)),
    ('mu', FieldEntry('mu', 'mu', 'cmb', (1, 1), 'OnGeoParamChanged', -1)),
    ('gamma', FieldEntry('gamma', 'gamma', 'txt', (0, 2), 'OnGeoParamChanged', -1)),
    ('b', FieldEntry('b', 'b', 'txt', (1, 2), 'OnGeoParamChanged', -1)),
    ('alpha', FieldEntry('alpha', 'alpha', 'txt', (0, 3), 'OnGeoParamChanged', -1)),
    ('d', FieldEntry('d', 'd', 'txt', (1, 3), 'OnGeoParamChanged', -1)),
    ('theta', FieldEntry('theta', 'theta', 'txt', (0, 4), 'OnGeoParamChanged', -1)),
    ('r', FieldEntry('r', 'r', 'txt', (1, 4), 'OnGeoParamChanged', -1))
])
# gravity component params
GRAVITY_CMPNTS = OrderedDict([
    ('gx', FieldEntry('GX', 'GX', 'txt', (0, 0), 'OnBaseTwistChanged', 0)),
    ('gy', FieldEntry('GY', 'GY', 'txt', (0, 1), 'OnBaseTwistChanged', 1)),
    ('gz', FieldEntry('GZ', 'GZ', 'txt', (0, 2), 'OnBaseTwistChanged', 2))
])
# robot type params
ROBOT_TYPE = OrderedDict([
    ('name', FieldEntry('Name of the robot:', 'name', 'lbl', (0, 0), None, -1)),
    ('num_links', FieldEntry('Number of moving links:', 'NL', 'lbl', (1, 0), None, -1)),
    ('num_joints', FieldEntry('Number of joints:', 'NJ', 'lbl', (2, 0), None, -1)),
    ('num_frames', FieldEntry('Number of frames:', 'NF', 'lbl', (3, 0), None, -1)),
    ('structure', FieldEntry('Type of structure:', 'type', 'lbl', (4, 0), None, -1)),
    ('is_floating', FieldEntry('Is Floating Base:', 'floating', 'lbl', (5, 0), None, -1)),
    ('is_mobile', FieldEntry('Is Mobile Robot:', 'mobile', 'lbl', (6, 0), None, -1)),
    ('num_loops', FieldEntry('Number of closed loops:', 'loops', 'lbl', (7, 0), None, -1))
])

# menu bar
# TODO: add event handlers to dict entries as well
MAIN_MENU = OrderedDict([
    ('file_menu', '&File'),
    ('geom_menu', '&Geometric'),
    ('kin_menu', '&Kinematic'),
    ('dyn_menu', '&Dynamic'),
    ('iden_menu', '&Identification'),
    ('optim_menu', '&Optimiser'),
    ('viz_menu', '&Visualisation'),
])
VIZ_MENU = OrderedDict(m_viz='&Visualisation')
IDEN_MENU = OrderedDict([
    ('m_base_inertial_params', 'Base Inertial parameters'),
    ('m_dyn_iden_model', 'Dynamic Identification Model'),
    ('m_energy_iden_model', 'Energy Identification Model'),
])
DYN_MENU = OrderedDict([
    ('m_idym', 'Inverse Dynamic Model'),
    ('m_inertia_matrix', 'Inertia matrix'),
    ('m_h_term', 'Centrifugal, Coriolis & Gravity torques'),
    ('m_ddym', 'Direct Dynamic Model'),
])
KIN_MENU = OrderedDict([
    ('m_jac_matrix', 'Jacobian matrix'),
    ('m_determinant', 'Determinant of a Jacobian'),
    ('m_kin_constraint', 'Kinematic constraint equation of loops'),
    ('m_vel', 'Velocities'),
    ('m_acc', 'Accelerations'),
    ('m_jpqp', 'Jpqp'),
])
GEOM_MENU = OrderedDict([
    ('m_trans_matrix', 'Transformation matrix'),
    ('m_fast_dgm', 'Fast Geometric model'),
    ('m_igm_paul', 'IGM - Paul method'),
    ('m_igm_pieper', 'IGM - Pieper method'),
    ('m_geom_constraint', 'Geometric constraint equation of loops'),
])
FILE_MENU = OrderedDict([
    ('m_new', '&New'),
    ('m_open', '&Open'),
    ('m_save', '&Save'),
    ('m_save_as', 'Save &As New Robot'),
    ('m_pref', 'Preferences -- (unavailable)'),
    ('m_exit', 'E&xit'),
])

