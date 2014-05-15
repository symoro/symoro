# -*- coding: utf-8 -*-


"""
This module contains the functions and classes for the computation of
dynamic models (inverse and direct).
"""


def _compute_link_velocity(robo, j, i):
    if i == 0: robo.vels.append(robo.base_vel.val)
    j_s_i = robo.geos[j].tmat.s_i_wrt_j
    i_v_i = robo.vels[i]
    qdot = robo.qdots[j]
    return (j_s_i * i_v_i) + (qdot * robo.geos[j].axisa)
    pass


def inverse_dynamic_model(robo):
    """
    Compute the inverse dynamic model for the given robot by using the
    recursive Newton-Euler algorithm.

    Args:
        robo: An instance of the FloatingRobot class.

    Returns:
        The inverse dynamic model of the robot.
    """
    # some new attributes for the robot
    robo.vels = list()
    # some book keeping variables
    nums = robo.joint_nums
    # first forward recursion
    for j in nums:
        if j == 0: continue
        # antecedent index
        i = robo.geos[j].ant
        # compute j^S_i : screw transformation matrix
        j_s_i = robo.geos[j].tmat.s_i_wrt_j
        # compute j^V_j : link velocity (6x1)
        # compute j^gamma_j : gyroscopic acceleration (6x1)
        # compute j^beta_j : external+coriolis+centrifugal wrench (6x1)
        # compute j^zeta_j : relative acceleration (6x1)
        pass
    # backward recursion
    for j in reversed(nums):
        if j == 0: continue
        pass
    # second forward recursion
    for j in robo.nums:
        if j == 0: continue
        pass
    pass

