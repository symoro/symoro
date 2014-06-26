# -*- coding: utf-8 -*-


"""
This module contains the methods used to initialise the different
matrices and parameters for various models.
"""

from copy import copy

from sympy import zeros
from sympy import Matrix
from symoroutils import tools


class ParamsInit(object):
    """
    This class contains methods that are used to initialise the different
    matrices and parameters for various models. All the methods in this
    class are class-methods.
    """
    @classmethod
    def init_jplus(cls, robo):
        """Copies the inertia parameters.
        Used for composed link inertia computation

        Returns
        =======
        Jplus: list of Matrices 3x3
        MSplus: list of Matrices 3x1
        Mplus: list of var
        """
        j_plus = copy(robo.J)
        j_plus.append(zeros((3, 3)))
        ms_plus = copy(robo.MS)
        ms_plus.append(zeros((3, 1)))
        m_plus = copy(robo.M)
        m_plus.append(0)
        return j_plus, ms_plus, m_plus

    @classmethod
    def init_mat(cls, robot, num=3):
        """Generates a list of zero Matrices.

        Generates a list of robot.NL Matrices of size num x num. robot.NL is
        the robot's number of links.

        Parameters
        ==========
        robot: Robot
            Instance of robot description container
        num: int, optional
            size of the matries, default is 3

        Returns
        =======
        list of Matrices num x num
        """
        return [zeros((num, num))] * (robot.NL)

    @classmethod
    def init_vec(cls, robot, num=3, ext=0):
        """Generates a list of vectors.

        The size of the list is the number of robot's links (robot.NL) by
        default. Actually, it is (robot.NL + ext). Each vector has size
        num x 1.

        Parameters
        ==========
        robot: Robot
            Instance of robot description container
        num: int, optional
            size of the vectors, default is 3
        ext: int, optional
            additional vector instances over number of links

        Returns
        =======
        list of Matrices Nx1
        """
        return [zeros((num, 1))] * (robot.NL + ext)

    @classmethod
    def init_scalar(cls, robot):
        """Generates a list of scalars.

        The size of the list is the number of robot's links.
        """
        return [0] * (robot.NL)

    @classmethod
    def init_w(cls, robot):
        """Generates a list of vectors for angular velocities.

        The size of the list is the number of robot's links.
        The zero vector is the base angular velocity
        """
        # TODO: the old help stated that list length should be robot.NL + 1
        omega = cls.init_vec(robot)
        omega[0] = robot.w0
        return omega

    @classmethod
    def init_v(cls, robot):
        """Generates a list of vectors for linear velocities.

        The size of the list is the number of robot's links.
        The zero vector is the base angular velocity.
        """
        # TODO: the old help stated that list length should be robot.NL + 1
        vel = cls.init_vec(robot)
        vel[0] = robot.v0
        return vel

    @classmethod
    def init_wv_dot(cls, robot, gravity=True):
        """Generates lists of vectors for angular and linear accelerations.

        The size of the list is the number of robot's links.
        The zero vector is the base angular velocity.

        Returns
        =======
        vdot: list of Matrices 3x1
        wdot: list of Matrices 3x1
        """
        # TODO: the old help stated that list length should be robot.NL + 1
        wdot = cls.init_vec(robot)
        wdot[0] = robot.wdot0
        vdot = cls.init_vec(robot)
        vdot[0] = robot.vdot0
        if gravity:
            vdot[0] -= robot.G
        return wdot, vdot

    @classmethod
    def init_u(cls, robot):
        """Generates a list of auxiliary U matrices"""
        u_aux_matrix = ParamsInit.init_mat(robot)
        # the value for the -1th base frame.
        u_aux_matrix.append(
                tools.skew(robot.w0)**2 + tools.skew(robot.wdot0)
            )
        return u_aux_matrix

    @classmethod
    def product_combinations(cls, vec):
        """Generates 6-vector of different v elements'
        product combinations

        Parameters
        ==========
        vec: Matrix 3x1
            vector

        Returns
        =======
        product_combinations: Matrix 6x1
        """
        return Matrix([
            vec[0]*vec[0], vec[0]*vec[1], vec[0]*vec[2],
            vec[1]*vec[1], vec[1]*vec[2], vec[2]*vec[2]
            ])
