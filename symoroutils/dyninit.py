# -*- coding: utf-8 -*-


"""
This module contains the parameters used to initialise the dynamic
model of a robot.
"""


class Init:
    @classmethod
    def init_Jplus(cls, robo):
        """Copies the inertia parameters.
        Used for composed link inertia computation

        Returns
        =======
        Jplus: list of Matrices 3x3
        MSplus: list of Matrices 3x1
        Mplus: list of var
        """
        Jplus = copy(robo.J)
        Jplus.append(zeros(3, 3))
        MSplus = copy(robo.MS)
        MSplus.append(zeros(3, 1))
        Mplus = copy(robo.M)
        Mplus.append(0)
        return Jplus, MSplus, Mplus

    @classmethod
    def init_mat(cls, robo, N=3):
        """Generates a list of Matrices.Size of the
        list is number of links.

        Parameters
        ==========
        robo: Robot
            Instance of robot description container
        N: int, optional
            size of the matries, default is 3

        Returns
        =======
        list of Matrices NxN
        """
        return [zeros(N, N) for i in xrange(robo.NL)]

    @classmethod
    def init_vec(cls, robo, N=3, ext=0):
        """Generates a list of vectors.
        Size of the list is number of links.

        Parameters
        ==========
        robo: Robot
            Instance of robot description container
        N: int, optional
            size of the vectors, default is 3
        ext: int, optional
            additional vector instances over number of links

        Returns
        =======
        list of Matrices Nx1
        """
        return [zeros(N, 1) for i in xrange(robo.NL+ext)]

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
        w = cls.init_vec(robo)
        w[0] = robo.w0
        return w

    @classmethod
    def init_v(cls, robo):
        """Generates a list of vectors for linear velocities.
        Size of the list is number of links + 1.
        The zero vector is the base angular velocity
        """
        v = cls.init_vec(robo)
        v[0] = robo.v0
        return v

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
    def init_U(cls, robo):
        """Generates a list of auxiliary U matrices"""
        U = Init.init_mat(robo)
        # the value for the -1th base frame
        U.append(tools.skew(robo.w0)**2 + tools.skew(robo.wdot0))
        return U

    @classmethod
    def product_combinations(cls, v):
        """Generates 6-vector of different v elements'
        product combinations

        Parameters
        ==========
        v: Matrix 3x1
            vector

        Returns
        =======
        product_combinations: Matrix 6x1
        """
        return Matrix([v[0]*v[0], v[0]*v[1], v[0]*v[2],
                       v[1]*v[1], v[1]*v[2], v[2]*v[2]])
