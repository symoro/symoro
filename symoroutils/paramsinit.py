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
        j_plus.append(zeros(3, 3))
        ms_plus = copy(robo.MS)
        ms_plus.append(zeros(3, 1))
        m_plus = copy(robo.M)
        m_plus.append(0)
        return j_plus, ms_plus, m_plus

    @classmethod
    def init_mat(cls, robo, num=3):
        """Generates a list of Matrices.Size of the
        list is number of links.

        Parameters
        ==========
        robo: Robot
            Instance of robot description container
        num: int, optional
            size of the matries, default is 3

        Returns
        =======
        list of Matrices numxnum
        """
        return [zeros(num, num) for i in xrange(robo.NL)]

    @classmethod
    def init_vec(cls, robo, num=3, ext=0):
        """Generates a list of vectors.
        Size of the list is number of links.

        Parameters
        ==========
        robo: Robot
            Instance of robot description container
        num: int, optional
            size of the vectors, default is 3
        ext: int, optional
            additional vector instances over number of links

        Returns
        =======
        list of Matrices Nx1
        """
        return [zeros(num, 1) for i in xrange(robo.NL+ext)]

    @classmethod
    def init_scalar(cls, robo):
        """Generates a list of vars.
        Size of the list is number of links.
        """
        return [0 for i in xrange(robo.NL)]

    @classmethod
    def init_w(cls, robo):
        """Generates a list of vectors for angular velocities.
        Size of the list is number of links + 1.
        The zero vector is the base angular velocity
        """
        omega = cls.init_vec(robo)
        omega[0] = robo.w0
        return omega

    @classmethod
    def init_v(cls, robo):
        """Generates a list of vectors for linear velocities.
        Size of the list is number of links + 1.
        The zero vector is the base angular velocity
        """
        vel = cls.init_vec(robo)
        vel[0] = robo.v0
        return vel

    @classmethod
    def init_wv_dot(cls, robo, gravity=True):
        """Generates lists of vectors for
        angular and linear accelerations.
        Size of the list is number of links + 1.
        The zero vector is the base angular velocity

        Returns
        =======
        vdot: list of Matrices 3x1
        wdot: list of Matrices 3x1
        """
        wdot = cls.init_vec(robo)
        wdot[0] = robo.wdot0
        vdot = cls.init_vec(robo)
        vdot[0] = robo.vdot0
        if gravity:
            vdot[0] -= robo.G
        return wdot, vdot

    @classmethod
    def init_u(cls, robo):
        """Generates a list of auxiliary U matrices"""
        u_aux_matrix = ParamsInit.init_mat(robo)
        # the value for the -1th base frame
        u_aux_matrix.append(
            tools.skew(robo.w0)**2 + tools.skew(robo.wdot0)
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


