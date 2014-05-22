# -*- coding: utf-8 -*-


"""
This module contains the functions and classes for the computation of
dynamic models (inverse and direct).
"""


from sympy import Matrix
from sympy import sign

from pysymoro.screw import Screw
from pysymoro.screw6 import Screw6
from symoroutils.tools import skew


class DynModel(object):
    """
    Data structure:
        Hold the various components obtained during the computation of
        dynamic models (inverse and direct).
    """
    def __init__(self, joints):
        """
        Constructor period.

        Args:
            joints: An iterative object containing all the joint
                numbers. This is usually the `joint_nums` attribute in
                `Robot` class.
        """
        # link velocity
        self.vels = list(None for j in joints)
        # gyroscopic acceleration
        self.gammas = list(None for j in joints)
        # wrench - external+coriolis+centrifugal
        self.betas = list(None for j in joints)
        # relative acceleration
        self.zetas = list(None for j in joints)
        # composite spatial inertia matrix
        self.composite_inertias = list(None for j in joints)
        # composite wrench
        self.composite_betas = list(None for j in joints)
        # link acceleration
        self.accels = list(None for j in joints)
        # reaction wrench
        self.wrenchs = list(None for j in joints)
        # joint torque
        self.torques = list(None for j in joints)

    def __str__(self):
        str_format = ""
        # add header
        str_format = str_format + "DynModel:\n"
        str_format = str_format + "---------\n"
        # get all the attributes currently in the class
        attrs = [attr for attr in dir(self) if not attr.startswith('_')]
        # add each attribute
        for attr in attrs:
            items = getattr(self, attr)
            if hasattr(items, '__iter__'):
                attr_str = self._str_items(items)
            else:
                attr_str = str(items) + '\n'
            str_format = str_format + str(attr) + ": \n"
            str_format = str_format + attr_str
        return str_format

    def __repr__(self):
        return str(self)

    def _str_items(self, items):
        """
        Create a string representation for a given attribute.

        Args:
            items: An attribute of the class which is a list

        Returns:
            A string representation of the attribute.
        """
        row_format = '\t' + ('{0:^6} : {1}') + '\n'
        str_format = ""
        for idx, item in enumerate(items):
            str_format = str_format + row_format.format(*(
                str(idx), str(item)
            ))
        return str_format


def _init_composite_inertia(model, robo, j):
    """
    Initialise the composite spatial inertia matrix of link j.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number

    Returns:
        An instance of DynModel that contains all the new values.
    """
    model.composite_inertias[j] = robo.dyns[j].spatial_inertia
    return model


def _init_composite_beta(model, robo, j):
    """
    Initialise the composite beta wrench of link j.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number

    Returns:
        An instance of DynModel that contains all the new values.
    """
    model.composite_betas[j] = model.betas[j]
    return model


def _compute_link_velocity(model, robo, j, i):
    """
    Compute the velocity of link j whose antecedent is i.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number
        i: antecendent value

    Returns:
        An instance of DynModel that contains all the new values.
    """
    j_v_j = Screw()
    if i == 0: model.vels[i] = robo.base_vel
    # local variables
    j_s_i = robo.geos[j].tmat.s_i_wrt_j
    qdot_j = robo.qdots[j]
    j_a_j = robo.geos[j].axisa
    i_v_i = model.vels[i].val
    # actual computation
    j_v_j.val = (j_s_i * i_v_i) + (qdot_j * j_a_j)
    # store computed velocity in model
    model.vels[j] = j_v_j
    return model


def _compute_link_acceleration(model, robo, j, i):
    """
    Compute the acceleration of link j whose antecedent is i.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number
        i: antecendent value

    Returns:
        An instance of DynModel that contains all the new values.
    """
    j_vdot_j = Screw()
    if i == 0 and not robo.is_floating:
        model.accels[i] = robo.base_accel
    # local variables
    j_s_i = robo.geos[j].tmat.s_j_wrt_i
    i_vdot_i = model.accels[i].val
    j_zeta_j = model.zetas[j].val
    # actual computation
    j_vdot_j.val = (j_s_i * i_vdot_i) + j_zeta_j
    # store computed velocity in model
    model.accels[j] = j_vdot_j
    return model


def _compute_gyroscopic_acceleration(model, robo, j, i):
    """
    Compute the gyroscopic acceleration of link j.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number
        i: antecendent value

    Returns:
        An instance of DynModel that contains all the new values.
    """
    j_gamma_j = Screw()
    # local variables
    j_rot_i = robo.geos[j].tmat.inv_rot
    i_trans_j = robo.geos[j].tmat.trans
    i_omega_i = model.vels[i].ang
    sigma_j = robo.geos[j].sigma
    sigma_dash_j = 1 - sigma_j
    j_z_j = robo.geos[j].zunit
    qdot_j = robo.qdots[j]
    # actual computation
    j_omega_i = j_rot_i * i_omega_i
    # term1 = i_omega_i x (i_omega_i x i_trans_j)
    term1 = skew(i_omega_i) * (skew(i_omega_i) * i_trans_j)
    # term2 = j_omega_i x (qdot_j * j_z_j)
    term2 = skew(j_omega_i) * (qdot_j * j_z_j)
    j_gamma_j.lin = (j_rot_i * term1) + (2 * sigma_j * term2)
    j_gamma_j.ang = sigma_dash_j * term2
    # store computed acceleration in model
    model.gammas[j] = j_gamma_j
    return model


def _compute_relative_acceleration(model, robo, j):
    """
    Compute the relative acceleration of link j.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number

    Returns:
        An instance of DynModel that contains all the new values.
    """
    j_zeta_j = Screw()
    # local variables
    j_a_j = robo.geos[j].axisa
    qddot_j = robo.qddots[j]
    j_gamma_j = model.gammas[j].val
    # actual computation
    j_zeta_j.val = j_gamma_j + (qddot_j * j_a_j)
    # store computed relative acceleration in model
    model.zetas[j] = j_zeta_j
    return model


def _compute_beta_wrench(model, robo, j):
    """
    Compute the wrench for link j which combines the external forces,
    Coriolis forces and centrifugal forces.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number

    Returns:
        An instance of DynModel that contains all the new values.
    """
    j_beta_j = Screw()
    # local variables
    j_omega_j = model.vels[j].ang
    j_fe_j = robo.dyns[j].wrench.val
    j_ms_j = robo.dyns[j].mass_tensor
    j_inertia_j = robo.dyns[j].inertia
    # actual computation
    # lin_term = j_omega_j x (j_omega_j x j_ms_j)
    lin_term = skew(j_omega_j) * (skew(j_omega_j) * j_ms_j)
    # ang_term = j_omega_j x (j_inertia_j * j_omega_j)
    ang_term = skew(j_omega_j) * (j_inertia_j * j_omega_j)
    term = Screw(lin=lin_term, ang=ang_term)
    j_beta_j.val = - j_fe_j - term.val
    # store computed wrench in model
    model.betas[j] = j_beta_j
    return model


def _compute_composite_inertia(model, robo, j, i):
    """
    Compute the composite spatial inertia matrix for link i.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number
        i: antecedent value

    Returns:
        An instance of DynModel that contains all the new values.
    """
    i_inertia_i_c = Screw6()
    # local variables
    j_s_i = robo.geos[j].tmat.s_i_wrt_j
    i_inertia_i = model.composite_inertias[i].val
    j_inertia_j_c = model.composite_inertias[j].val
    # actual computation
    i_inertia_i_c.val = i_inertia_i + \
        (j_s_i.transpose() * j_inertia_j_c * j_s_i)
    # store computed matrix in model
    model.composite_inertias[i] = i_inertia_i_c
    return model


def _compute_composite_beta(model, robo, j, i):
    """
    Compute the composite beta wrench for link i.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number
        i: antecedent value

    Returns:
        An instance of DynModel that contains all the new values.
    """
    i_beta_i_c = Screw()
    # local variables
    j_s_i = robo.geos[j].tmat.s_i_wrt_j
    i_beta_i = model.composite_betas[i].val
    j_beta_j_c = model.composite_betas[j].val
    j_inertia_j_c = model.composite_inertias[j].val
    j_zeta_j = model.zetas[j].val
    # actual computation
    i_beta_i_c.val = i_beta_i - (j_s_i.transpose() * j_beta_j_c) + \
        (j_s_i.transpose() * j_inertia_j_c * j_zeta_j)
    # store computed beta in model
    model.composite_betas[i] = i_beta_i_c
    return model


def _compute_reaction_wrench(model, robo, j):
    """
    Compute the reaction wrench for link j.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: link number

    Returns:
        An instance of DynModel that contains all the new values.
    """
    j_f_j = Screw()
    # local variables
    j_vdot_j = model.accels[j].val
    j_inertia_j_c = model.composite_inertias[j].val
    j_beta_j_c = model.composite_betas[j].val
    # actual computation
    j_f_j.val = (j_inertia_j_c * j_vdot_j) - j_beta_j_c
    # store computed reaction wrench in model
    model.wrenchs[j] = j_f_j
    return model


def _compute_joint_torque(model, robo, j):
    """
    Compute the joint torque for joint j.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot
        j: joint number

    Returns:
        An instance of DynModel that contains all the new values.
    """
    # local variables
    qdot_j = robo.qdots[j]
    qddot_j = robo.qddots[j]
    j_a_j = robo.geos[j].axisa
    ia_j = robo.dyns[j].ia
    f_cj = robo.dyns[j].frc
    f_vj = robo.dyns[j].frv
    j_f_j = model.wrenchs[j].val
    # actual computation
    wrench_term = j_f_j.transpose() * j_a_j
    actuator_inertia_term = Matrix([ia_j * qddot_j])
    coriolis_friction_term = Matrix([f_cj * sign(qdot_j)])
    viscous_friction_term = Matrix([f_vj * qdot_j])
    gamma_j = wrench_term + actuator_inertia_term + \
        viscous_friction_term + coriolis_friction_term
    # store computed torque in model
    model.torques[j] = gamma_j
    return model


def _compute_base_acceleration(model, robo):
    """
    Compute the base acceleration for a robot with floating base without
    and with taking gravity into account. In the case of a robot with
    fixed base, this function returns just the effect of gravity as the
    base acceleration.

    Args:
        model: An instance of DynModel
        robo: An instance of Robot

    Returns:
        An instance of DynModel that contains all the new values.
    """
    o_vdot_o = Screw()
    gravity = Screw()
    # local variables
    gravity.lin = robo.gravity
    o_inertia_o_c = model.composite_inertias[0].val
    o_beta_o_c = model.composite_betas[0].val
    # actual computation
    if robo.is_floating:
        # TODO: replace sympy's matrix inversion with custom function
        o_vdot_o.val = o_inertia_o_c.inv() * o_beta_o_c
    # store computed base acceleration without gravity effect in model
    model.base_accel_no_gravity = o_vdot_o
    # compute base acceleration taking gravity into account
    o_vdot_o.val = o_vdot_o.val + gravity.val
    # store in model
    model.accels[0] = o_vdot_o
    return model


def inverse_dynamic_model(robo):
    """
    Compute the inverse dynamic model for the given robot by using the
    recursive Newton-Euler algorithm.

    Args:
        robo: An instance of the FloatingRobot class.

    Returns:
        The inverse dynamic model of the robot.
    """
    # some book keeping variables
    model = DynModel(robo.joint_nums)
    # first forward recursion
    for j in robo.joint_nums:
        if j == 0: continue
        # antecedent index
        i = robo.geos[j].ant
        # compute j^V_j : link velocity (6x1)
        model = _compute_link_velocity(model, robo, j, i)
        # compute j^gamma_j : gyroscopic acceleration (6x1)
        model = _compute_gyroscopic_acceleration(model, robo, j, i)
        # compute j^beta_j : external+coriolis+centrifugal wrench (6x1)
        model = _compute_beta_wrench(model, robo, j)
        # compute j^zeta_j : relative acceleration (6x1)
        # TODO: check joint flexibility
        model = _compute_relative_acceleration(model, robo, j)
    # first backward recursion - initialisation step
    for j in reversed(robo.joint_nums):
        if j == 0:
            # compute 0^beta_0
            model = _compute_beta_wrench(model, robo, j)
        # initialise j^I_j^c : composite spatial inertia matrix
        model = _init_composite_inertia(model, robo, j)
        # initialise j^beta_j^c : composite wrench
        model = _init_composite_beta(model, robo, j)
    # second backward recursion - compute composite terms
    for j in reversed(robo.joint_nums):
        if j == 0:
            # compute 0^\dot{V}_0 : base acceleration
            # for fixed base robots, the value returned is just the
            # effect of gravity
            model = _compute_base_acceleration(model, robo)
            continue
        # antecedent index
        i = robo.geos[j].ant
        # compute i^I_i^c : composite spatial inertia matrix
        model = _compute_composite_inertia(model, robo, j, i)
        # compute i^beta_i^c : composite wrench
        model = _compute_composite_beta(model, robo, j, i)
    # second forward recursion
    for j in robo.joint_nums:
        if j == 0: continue
        # antecedent index
        i = robo.geos[j].ant
        # compute j^\dot{V}_j : link acceleration
        model = _compute_link_acceleration(model, robo, j, i)
        # compute j^F_j : reaction wrench
        model = _compute_reaction_wrench(model, robo, j)
        # compute gamma_j : joint torque
        model = _compute_joint_torque(model, robo, j)
    return model


def direct_dynamic_model(robo):
    """
    Compute the direct dynamic model for the given robot by using the
    recursive Newton-Euler algorithm.

    Args:
        robo: An instance of the FloatingRobot class.

    Returns:
        The direct dynamic model of the robot.
    """
    # some book keeping variables
    model = DynModel(robo.joint_nums)
    # TODO: complete direct dynamic model
    return model


