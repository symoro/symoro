# -*- coding: utf-8 -*-


"""
This module contains the functions and classes for the computation of
dynamic models (inverse and direct).
"""


from pysymoro.screw import Screw
from symoroutil.tools import skew


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
        self.vels = list(None for j in joints)
        self.gammas = list(None for j in joints)
        self.betas = list(None for j in joints)
        self.zetas = list(None for j in joints)

    def __str__(self):
        str_format = ""
        # add header
        str_format = str_format + "DynModel:\n"
        str_format = str_format + "---------\n"
        # add link velocity
        str_format = str_format + "vels: \n"
        str_format = str_format + self._str_items(self.vels)
        # add gyroscopic acceleration
        str_format = str_format + "gammas: \n"
        str_format = str_format + self._str_items(self.gammas)
        # add wrench - external + coriolis + centrifugal
        str_format = str_format + "betas: \n"
        str_format = str_format + self._str_items(self.betas)
        # add relative acceleration
        str_format = str_format + "zetas: \n"
        str_format = str_format + self._str_items(self.zetas)
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
    i_v_i = models.vels[i].val
    qdot_j = robo.qdots[j]
    j_a_j = robo.geos[j].axisa
    # actual computation
    j_v_j.val = (j_s_i * i_v_i) + (qdot_j * j_a_j)
    # store computed velocity in model
    model.vels[j] = j_v_j
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
    j_rot_i = robo.geos[j].tmat.rot_inv
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
    nums = robo.joint_nums
    model = DynModel(nums)
    # first forward recursion
    for j in nums:
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
    # backward recursion
    for j in reversed(nums):
        # antecedent index
        i = robo.geos[j].ant
        if j != 0:
            pass
        else:
            pass
    # second forward recursion
    for j in robo.nums:
        if j == 0: continue
        pass
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
    pass


