# -*- coding: utf-8 -*-

"""
This module of SYMORO package provides description
of the robot parametrizaion container and symbol replacer class.

The core symbolic library is sympy.

ECN - ARIA1 2013
"""


import re
import os
from copy import copy
from itertools import combinations

from sympy import sin, cos, sign, pi
from sympy import Symbol, Matrix, Expr, Integer
from sympy import Mul, Add, factor, zeros, var, sympify, eye

from symoroutils import filemgr


ZERO = Integer(0)
ONE = Integer(1)
CLOSED_LOOP = 'Closed loop'
SIMPLE = 'Simple'
TREE = 'Tree'
TYPES = [SIMPLE, TREE, CLOSED_LOOP]
FAIL = 1
OK = 0
INT_KEYS = ['ant', 'sigma', 'mu']


#TODO: write consistency check
#TODO: Ask about QP QDP file writing. Number of joints is different
#from number of links
class Robot:
    """Container of the robot parametric description.
    Responsible for low-level geometric transformation
    and direct geometric model generation.
    Also provides different representations of parameters."""
    def __init__(self, name, NL=0, NJ=0, NF=0, is_mobile=False,
                 structure=TREE):
        # member variables:
        self.name = name
        """  name of the robot: string"""
        self.directory = filemgr.get_folder_path(name)
        """ directory name"""
        self.is_mobile = is_mobile
        """ whethere the base frame is floating: bool"""
        self.nl = NL
        """  number of links: int"""
        self.nj = NJ
        """  number of joints: int"""
        self.nf = NF
        """  number of frames: int"""
        self.structure = structure
        """ type of robot's structure"""
        self.sigma = [0 for i in xrange(NF + 1)]
        """  joint type: list of int"""
        self.ant = range(-1, self.NF - 1)
        """  index of antecedent joint: list of int"""
        self.mu = [0 for i in xrange(NF + 1)]
        """motorization, if 1, then the joint im motorized"""
        self.theta = [0] + [var('th%s' % (i+1)) for i in xrange(NF)]
        """  geometrical parameter: list of var"""
        self.r = [0 for i in xrange(NF + 1)]
        """  geometrical parameter: list of var"""
        self.alpha = [0 for i in xrange(NF + 1)]
        """  geometrical parameter: list of var"""
        self.d = [0 for i in xrange(NF + 1)]
        """  geometrical parameter: list of var"""
        self.gamma = [0 for i in xrange(NF + 1)]
        """  geometrical parameter: list of var"""
        self.b = [0 for i in xrange(NF + 1)]
        """  geometrical parameter: list of var"""
        self.Z = eye(4)
        """ transformation from reference frame to zero frame"""
        num = range(self.NL)
        numj = range(self.NJ)
        self.w0 = zeros(3, 1)
        """  base angular velocity: 3x1 matrix"""
        self.wdot0 = zeros(3, 1)
        """  base angular acceleration: 3x1 matrix"""
        self.v0 = zeros(3, 1)
        """  base linear velocity: 3x1 matrix"""
        self.vdot0 = zeros(3, 1)
        """  base linear acceleration: 3x1 matrix"""
        self.qdot = [var('QP{0}'.format(i)) for i in numj]
        """  joint speed: list of var"""
        self.qddot = [var('QDP{0}'.format(i)) for i in numj]
        """  joint acceleration: list of var"""
        self.Nex = [zeros(3, 1) for i in num]
        """  external moment of link: list of 3x1 matrix"""
        self.Nex[-1] = Matrix(var('CX{0}, CY{0}, CZ{0}'.format(self.NL - 1)))
        self.Fex = [zeros(3, 1) for i in num]
        """  external force of link: list of 3x1 matrix"""
        self.Fex[-1] = Matrix(var('FX{0}, FY{0}, FZ{0}'.format(self.NL - 1)))
        self.FS = [var('FS{0}'.format(i)) for i in num]
        """  dry friction coefficient: list of ver"""
        self.IA = [var('IA{0}'.format(i)) for i in num]
        """  joint actuator inertia: list of var"""
        self.FV = [var('FV{0}'.format(i)) for i in num]
        """  viscous friction coefficient: list of var"""
        self.MS = [Matrix(var('MX{0}, MY{0}, MZ{0}'.format(i))) for i in num]
        """  first momentum of link: list of 3x1 matrix"""
        self.M = [var('M{0}'.format(i)) for i in num]
        """  mass of link: list of var"""
        self.GAM = [var('GAM{0}'.format(i)) for i in numj]
        """  joint torques: list of var"""
        J_str = 'XX{0},XY{0},XZ{0},XY{0},YY{0},YZ{0},XZ{0},YZ{0},ZZ{0}'
        self.J = [Matrix(3, 3, var(J_str.format(i))) for i in num]
        """  inertia tensor of link: list of 3x3 matrix"""
        self.G = Matrix([0, 0, var('G3')])
        """  gravity vector: 3x1 matrix"""

    # member methods:
    def put_val(self, j, name, val):
        #TODO: write proper parser
        #accepts tuple
        try:
            if isinstance(val, str) or isinstance(val, unicode):
                val = sympify(val)
                assert isinstance(val, Expr)
            geom_head = self.get_geom_head()
            if name in INT_KEYS:
                val = int(val)
        except:
            return FAIL
        base_vel_head = self.get_base_vel_head()
        ext_dynam_head = self.get_ext_dynam_head()
        dynam_head = self.get_dynam_head()
        ext_head = ext_dynam_head[7:] + ['IA']
        f_ex_head = ext_dynam_head[1:4]
        n_ex_head = ext_dynam_head[4:7]
        if name in ext_head + geom_head + base_vel_head:
            X = getattr(self, name)
            X[j] = val
        elif name in f_ex_head:
            self.Fex[j][f_ex_head.index(name)] = val
        elif name in n_ex_head:
            self.Nex[j][n_ex_head.index(name)] = val
        elif name in dynam_head:
            params = self.get_inert_param(j)
            i = dynam_head.index(name)
            params[i-1] = val
            self.put_inert_param(params, j)
        elif name == 'Z':
            self.Z[j] = val
        return OK

    def get_val(self, j, name):
        geom_head = self.get_geom_head()
        base_vel_head = self.get_base_vel_head()
        ext_dynam_head = self.get_ext_dynam_head()
        dynam_head = self.get_dynam_head()
        ext_head = ext_dynam_head[7:] + ['IA']
        f_ex_head = ext_dynam_head[1:4]
        n_ex_head = ext_dynam_head[4:7]
        if name in ext_head + geom_head + base_vel_head:
            X = getattr(self, name)
            return X[j]
        elif name in f_ex_head:
            return self.Fex[j][f_ex_head.index(name)]
        elif name in n_ex_head:
            return self.Nex[j][n_ex_head.index(name)]
        elif name in self.get_dynam_head():
            params = self.get_inert_param(j)
            i = dynam_head.index(name)
            return params[i-1]
        elif name == 'Z':
            return self.Z[j]

    def get_q_chain(self, j, k=0):
        """Generates vector of joint variables in chain
        between j, k (zero by default)
        """
        chain = self.chain(j, k)
        q = []
        for i in reversed(chain):
            if int(self.sigma[i]) == 0:
                q.append(self.theta[i])
            elif int(self.sigma[i]) == 1:
                q.append(self.r[i])
        return q

    def get_q(self, i):
        """ Returns symbol of joint variable
        or 0 if joint is fixed
        """
        if self.sigma[i] == 0:
            return self.theta[i]
        elif self.sigma[i] == 1:
            return self.r[i]
        else:
            return 0

    @property
    def q_vec(self):
        """Generates vector of joint variables
        """
        qs = []
        for i in xrange(1, self.NJ):
            if self.sigma[i] != 2:
                qs.append(self.get_q(i))
        return qs

    @property
    def endeffectors(self):
        return set(range(1, self.NJ + 1)) - set(self.ant)

    @property
    def q_passive(self):
        """Generates vector of passive joint variables (including cut!)
        """
        q = list()
        for i in xrange(1, self.NJ):
            if self.mu[i] == 0:
                q.append(self.get_q(i))
        return q

    @property
    def q_active(self):
        """Generates vector of active joint variables (including cut!)
        """
        q = list()
        for i in xrange(1, self.NJ):
            if self.mu[i] == 1:
                q.append(self.get_q(i))
        return q

    @property
    def indx_passive(self):
        """Generates vector of passive joint indices
        """
        return [i for i in xrange(1, self.NL) if self.mu[i] == 0]

    @property
    def indx_active(self):
        """Generates vector of active joint indices
        """
        return [i for i in xrange(1, self.NL) if self.mu[i] == 1]

    @property
    def indx_cut(self):
        """Generates vector of cut joint indices
        """
        return range(self.NL, self.NJ)

    def fric_v(self, j):
        """Fluid friction torque

        Parameters
        ==========
        j: int
            Joint index.

        Returns
        =======
        fric_v: sympy expression
            Expression for fluid friction torque of joint j
        """
        return self.FV[j] * self.qdot[j]

    def fric_s(self, j):
        """Dry friction torque

        Parameters
        ==========
        j: int
            Joint index.

        Returns
        =======
        fric_s: sympy expression
            Expression for dry friction torque of joint j
        """
        return self.FS[j] * sign(self.qdot[j])

    @property
    def W0(self):
        return self.w0

    @property
    def WP0(self):
        return self.wdot0

    @property
    def V0(self):
        return self.w0

    @property
    def VP0(self):
        return self.wdot0

    @property
    def QP(self):
        return self.qdot

    @property
    def QDP(self):
        return self.qddot

    @property
    def NJ(self):
        """ Actual number of joints counting 0
        """
        return self.nj + 1

    @property
    def NL(self):
        """ Actual number of links counting 0
        """
        return self.nl + 1

    @property
    def NF(self):
        """ Actual number of frames counting 0
        """
        return self.nf + 1

    @property
    def loop_terminals(self):
        B = self.NJ - self.NL
        return [(i, i+B) for i in xrange(self.NL, self.NJ)]

    def paral(self, i, j):
        if j is None:
            return False
        elif self.ant[i] == j:
            return sin(self.alpha[i]) == 0
        elif self.ant[j] == i:
            return sin(self.alpha[j]) == 0
        elif self.ant[j] == self.ant[i]:
            return sin(self.alpha[j] - self.alpha[i]) == 0
        else:
            return False

    def tau_ia(self, j):
        """Actuator inertia torque

        Parameters
        ==========
        j: int
            Joint index.

        Returns
        =======
        fric_v: sympy expression
            Expression for actuator inertia torque of joint j
        """
        return self.IA[j] * self.qddot[j]

    def get_angles(self, j):
        """List of non-constant angles of frame j

        Parameters
        ==========
        j: int
            Frame index.

        Returns
        =======
        get_angles: list of touples (var, name)
            Returns list of touples, where:
            var - the angle symbol,
            name - brief name for cos and sin abbreviation
        """
        angs = []
        if j not in xrange(self.NF):
            return angs
        if type(self.theta[j]) != int and not self.theta[j].is_number:
            angs.append((self.theta[j], j))
        if type(self.alpha[j]) != int and not self.alpha[j].is_number:
            angs.append((self.alpha[j], 'A%s' % j))
        if type(self.gamma[j]) != int and not self.gamma[j].is_number:
            angs.append((self.gamma[j], 'G%s' % j))
        return angs

    def chain(self, j, k=0):
        """Chain of antecedent frames between j-th and k-th frames

        Parameters
        ==========
        j: int
            Start frame index.
        k: int
            Final frame index.

        Returns
        =======
        u: list of ints
            List of antecedent frames. j is the first index in the list.
            k is not included
        """
        u = []
        while j != k and j != 0:
            u.append(j)
            j = self.ant[j]
        return u

    def loop_chain(self, i, j):
        k = self.common_root(i, j)
        chain = self.chain(i, k)
        chain.append(k)
        if k != j:
            chain.extend(reversed(self.chain(j, k)))
        return chain

    def common_root(self, i, j):
        """Common root j-th and i-th frames

        Parameters
        ==========
        j: int
            Frame index.
        i: int
            Frame index.

        Returns
        =======
        common_root: int
            The highest index of the common frame in chains for i and j.
            If they don't have common root, -1
        """
        u = self.chain(i)
        while True:
            if j in u or j == 0:
                return j
            j = self.ant[j]

    def get_inert_param(self, j):
        """Returns 10-vector of inertia paremeters of link j.

        Parameters
        ==========
        j: int
            Link index.

        Returns
        =======
        get_dynam_param: Matrix 10x1
        """
        K = [self.J[j][0], self.J[j][1], self.J[j][2], self.J[j][4],
             self.J[j][5], self.J[j][8], self.MS[j][0], self.MS[j][1],
             self.MS[j][2], self.M[j]]
        return Matrix(K)

    def put_inert_param(self, K, j):
        """Write the inertia parameters of link j from 10-vector K.

        Parameters
        ==========
        K: Matrix 10x1
            Vector of inertia parameters
        j: int
            Link index.
        """
        self.J[j] = Matrix([[K[0], K[1], K[2]],
                            [K[1], K[3], K[4]],
                            [K[2], K[4], K[5]]])
        self.MS[j] = Matrix(3, 1, K[6:9])
        self.M[j] = K[9]

    def get_ext_dynam_head(self):
        """Returns header for external forces and torques,
        friction parameters and joint speeds, accelerations.
        Used for output generation.

        Returns
        =======
        get_ext_dynam_head: list of strings
        """
        return ['j', 'FX', 'FY', 'FZ', 'CX', 'CY', 'CZ',
                'FS', 'FV', 'QP', 'QDP', 'GAM']

    def get_dynam_head(self):
        """Returns header for inertia parameters.
        Used for output generation.

        Returns
        =======
        get_dynam_head: list of strings
        """
        return ['j', 'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ',
                'MX', 'MY', 'MZ', 'M', 'IA']

    def get_geom_head(self):
        """Returns header for geometric parameters.
        Used for output generation.

        Returns
        =======
        get_geom_head: list of strings
        """
        return ['j', 'ant', 'sigma', 'mu', 'gamma', 'b',
                'alpha', 'd', 'theta', 'r']

    def get_base_vel_head(self):
        """Returns header for base velocities and gravity vector.
        Used for output generation.

        Returns
        =======
        get_base_vel_head: list of strings
        """
        return ['axis', 'W0', 'WP0', 'V0', 'VP0', 'G']

    def get_param_vec(self, head, j):
        params = list()
        axis_dict = {0: 'X', 1: 'Y', 2: 'Z'}
        for h in head:
            if h == 'j':
                params.append(j)
            elif h == 'axis':
                params.append(axis_dict[j])
            else:
                params.append(self.get_val(j, h))
        return params

    @classmethod
    def CartPole(cls):
        """Generates Robot instance of classical
        CartPole dynamic system.
        """
        #TODO: bring it to the new notation with 0-frame
        robo = Robot()
        robo.name = 'CartPole'
        robo.ant = (-1, 0)
        robo.sigma = (1, 0)
        robo.alpha = (pi/2, pi/2)
        robo.d = (0, 0)
        robo.theta = (pi/2, var('Th2'))
        robo.r = (var('R1'), 0)
        robo.b = (0, 0)
        robo.gamma = (0, 0)
        robo.num = range(1, 3)
        robo.NJ = 2
        robo.NL = 2
        robo.NF = 2
        robo.Nex = [zeros(3, 1) for i in robo.num]
        robo.Fex = [zeros(3, 1) for i in robo.num]
        robo.FS = [0 for i in robo.num]
        robo.IA = [0 for i in robo.num]
        robo.FV = [var('FV{0}'.format(i)) for i in robo.num]
        robo.MS = [zeros(3, 1) for i in robo.num]
        robo.MS[1][0] = var('MX2')
        robo.M = [var('M{0}'.format(i)) for i in robo.num]
        robo.GAM = [var('GAM{0}'.format(i)) for i in robo.num]
        robo.J = [zeros(3) for i in robo.num]
        robo.J[1][2, 2] = var('ZZ2')
        robo.G = Matrix([0, 0, -var('G3')])
        robo.w0 = zeros(3, 1)
        robo.wdot0 = zeros(3, 1)
        robo.v0 = zeros(3, 1)
        robo.vdot0 = zeros(3, 1)
        robo.q = var('R1, Th2')
        robo.qdot = var('R1d, Th2d')
        robo.qddot = var('R1dd, Th2dd')
        robo.num.append(0)
        return robo

    @classmethod
    def SR400(cls):
        #TODO: bring it to the new notation with 0-frame
        """Generates Robot instance of SR400"""
        robo = Robot('SR400', 8, 9, 10, False)
        robo.ant = [-1, 0, 1, 2, 3, 4, 5, 1, 7, 8, 3]
        robo.sigma = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2]
        robo.mu = [0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
        robo.alpha = [0, 0, -pi/2, 0, -pi/2, pi/2, -pi/2, -pi/2, 0, 0, 0]
        d_var = var('D:9')
        robo.d = [0, 0, d_var[2], d_var[3], d_var[4], 0, 0,
                  d_var[2], d_var[8], d_var[3], -d_var[8]]
        robo.theta = [0] + list(var('th1:10')) + [0]
        robo.r = [0, 0, 0, 0, var('RL4'), 0, 0, 0, 0, 0, 0]
        robo.b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        robo.gamma = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pi/2]
        robo.structure = CLOSED_LOOP
        return robo

    @classmethod
    def RX90(cls):
        """Generates Robot instance of RX90"""
        robo = Robot('RX90', 6, 6, 6, False)
        # table of geometric parameters RX90
        robo.sigma = [2, 0, 0, 0, 0, 0, 0, 0]
        robo.alpha = [0, 0, pi/2, 0, -pi/2, pi/2, -pi/2]
        robo.d = [0, 0, 0, var('D3'), 0, 0, 0]
        robo.theta = [0] + list(var('th1:7'))
        robo.r = [0, 0, 0, 0, var('RL4'), 0, 0]
        robo.b = [0, 0, 0, 0, 0, 0, 0]
        robo.gamma = [0, 0, 0, 0, 0, 0, 0]
        robo.mu = [0, 1, 1, 1, 1, 1, 1]
        robo.structure = SIMPLE
#        robo.w0 = zeros(3, 1)
#        robo.wdot0 = zeros(3, 1)
#        robo.v0 = zeros(3, 1)
#        robo.vdot0 = zeros(3, 1)
#        robo.qdot = [var('QP{0}'.format(i)) for i in num]
#        robo.qddot = [var('QDP{0}'.format(i)) for i in num]
#        robo.Nex= [zeros(3, 1) for i in num]
#        robo.Nex[-1] = Matrix(var('CX{0}, CY{0}, CZ{0}'.format(robo.NJ)))
#        robo.Fex = [zeros(3, 1) for i in num]
#        robo.Fex[-1] = Matrix(var('FX{0}, FY{0}, FZ{0}'.format(robo.NJ)))
#        robo.FS = [var('FS{0}'.format(i)) for i in num]
#        robo.IA = [var('IA{0}'.format(i)) for i in num]
#        robo.FV = [var('FV{0}'.format(i)) for i in num]
#        robo.MS = [Matrix(var('MX{0}, MY{0}, MZ{0}'.format(i))) for i in num]
#        robo.M = [var('M{0}'.format(i)) for i in num]
#        robo.GAM = [var('GAM{0}'.format(i)) for i in num]
#        robo.J = [Matrix(3, 3, var(('XX{0}, XY{0}, XZ{0}, '
#                            'XY{0}, YY{0}, YZ{0}, '
#                            'XZ{0}, YZ{0}, ZZ{0}').format(i))) for i in num]
#        robo.G = Matrix([0, 0, var('G3')])
#        robo.num.append(0)
        return robo


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
        U.append(hat(robo.w0)**2 + hat(robo.wdot0))
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


def hat(v):
    """skew-symmetry : Generates vectorial preproduct matrix

    Parameters
    ==========
    v: Matrix 3x1
        vector

    Returns
    =======
    hat: Matrix 3x3
    """
    return Matrix([[0, -v[2], v[1]],
                   [v[2], 0, -v[0]],
                   [-v[1], v[0], 0]])


def l2str(list_var, spacing=8):
    """Converts a list into string, that will be
    written into the text table.

    Parameters
    ==========
    list_var: list
        List to be converted
    spacing: int, optional
        Defines the size of one cell of the table

    Returns
    =======
    s: string
        String representation

    Notes
    =====
    l2str([1, 2, 3]) will be converted into '1      2      3      '
    """
    s = ''
    for i in list_var:
        s += str(i) + ' '*(spacing-len(str(i)))
    return s


def get_trig_couple_names(sym):
    names_s = find_trig_names(sym, r'S', 1)
    names_c = find_trig_names(sym, r'C', 1)
    return names_c & names_s


def find_trig_names(sym, pref=r'', pref_len=0, post=r'', post_len=0):
    search_res = re.findall(pref + r'[AGm0-9]*' + post, str(sym))
    if post_len == 0:
        return set([s[pref_len:] for s in search_res])
    else:
        return set([s[pref_len:-post_len] for s in search_res])


def get_max_coef_list(sym, x):
    return [get_max_coef_mul(s, x) for s in Add.make_args(sym)]


def get_max_coef(sym, x):
    return Add.fromiter(get_max_coef_mul(s, x) for s in Add.make_args(sym))


def get_max_coef_mul(sym, x):
    """
    """
    k, ex = x.as_coeff_Mul()
    coef = sym / k
    pow_x = ex.as_powers_dict()
    pow_c = coef.as_powers_dict()
    pow_c[-1] = 0
    for a, pa in pow_x.iteritems():
        na = -a
        if a in pow_c and pow_c[a] >= pa:
            pow_c[a] -= pa
        elif na in pow_c and pow_c[na] >= pa:
            pow_c[na] -= pa
            if pa % 2:
                pow_c[-1] += 1
        else:
            return ZERO
    return Mul.fromiter(c**p for c, p in pow_c.iteritems())


def ang_sum(np1, np2, nm1, nm2):
    np2, nm1 = reduce_str(np2, nm1)
    np1, nm2 = reduce_str(np1, nm2)
    if len(nm1) + len(nm2) == 0:
        return np1 + np2
    else:
        return np1 + np2 + 'm' + nm1 + nm2


def get_pos_neg(s):
    if s.find('m') != -1:
        s_split = s.split('m')
        return s_split[0], s_split[1]
    else:
        return s, ''


def reduce_str(s1, s2):
    while True:
        for j, char in enumerate(s1):
            if char in 'AG':
                i = s2.find(s1[j:j+2])
                k = 2
            else:
                i = s2.find(char)
                k = 1
            if i != -1:
                if i+k < len(s2):
                    s2_tail = s2[i+k:]
                else:
                    s2_tail = ''
                if j+k < len(s1):
                    s1_tail = s1[j+k:]
                else:
                    s1_tail = ''
                s2 = s2[:i] + s2_tail
                s1 = s1[:j] + s1_tail
                break
        else:
            break
    return s1, s2


def CS_syms(name):
    if isinstance(name, str) and name[0] == 'm':
        C, S = var('C{0}, S{0}'.format(name[1:]))
        return C, -S
    else:
        return var('C{0}, S{0}'.format(name))


def sym_less(A, B):
    A_measure = A.count_ops()
    B_measure = B.count_ops()
    return A_measure < B_measure


def get_angles(expr):
    angles_s = set()
    for s in expr.atoms(sin):
        angles_s |= set(s.args)
    angles_c = set()
    for c in expr.atoms(cos):
        angles_c |= set(c.args)
    return angles_s & angles_c


def cancel_terms(sym, X, coef):
    if coef.is_Add:
        for arg_c in coef.args:
            sym = cancel_terms(sym, X, arg_c)
    else:
        terms = Add.make_args(sym)
        return Add.fromiter(t for t in terms if t != X*coef)


def trigonometric_info(sym):
    if not sym.has(sin) and not sym.has(cos):
        short_form = True
        names = get_trig_couple_names(sym)
    else:
        short_form = False
        names = get_angles(sym)
    return names, short_form


