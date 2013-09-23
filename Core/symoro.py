"""
This module of SYMORO package provides description
of the robot parametrizaion container and symbol replacer class.

The core symbolic library is sympy.

ECN - ARIA1 2013
"""
import os
import re
from copy import copy
from itertools import combinations
from sympy import sin, cos, sign, pi
from sympy import Symbol, Matrix, Expr, Integer
from sympy import Mul, Add, factor, zeros, var, sympify, eye

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
        self.directory = name
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
        return [self.get_q(i) for i in xrange(1, self.NF)]

    @property
    def endeffectors(self):
        return set(range(1, self.NJ + 1)) - set(self.ant)

    @property
    def q_passive(self):
        """Generates vector of passive joint variables
        """
        q = list()
        for i in xrange(1, self.NJ):
            if self.mu[i] == 0:
                q.append(self.get_q(i))
        return q

    @property
    def q_active(self):
        """Generates vector of active joint variables
        """
        q = list()
        for i in xrange(1, self.NJ):
            if self.mu[i] == 1:
                q.append(self.get_q(i))
        return q

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
        The last vector is the base angular velocity
        """
        w = cls.init_vec(robo)
        w.append(robo.w0)
        return w

    @classmethod
    def init_wv_dot(cls, robo):
        """Generates lists of vectors for
        angular and linear accelerations.
        Size of the list is number of links + 1.
        The last vector is the base angular velocity

        Returns
        =======
        vdot: list of Matrices 3x1
        wdot: list of Matrices 3x1
        """
        wdot = cls.init_vec(robo)
        wdot.append(robo.wdot0)
        vdot = cls.init_vec(robo)
        vdot.append(robo.vdot0 - robo.G)
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
    """Generates vectorial preproduct matrix

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


def make_fname(robo, ext=None):
    if ext is None:
        fname = '%s.par' % robo.name
    else:
        fname = '%s_%s.txt' % (robo.name, ext)
    if robo.directory.find(':') == -1:      # path is not absolute
        full_name = 'Robots\\%s\\%s' % (robo.directory, fname)
        if not os.path.exists('Robots'):
            os.makedirs('Robots')
        d = 'Robots\\%s' % robo.directory
        if not os.path.exists(d):
            os.makedirs(d)
    else:
        full_name = '%s\\%s' % (robo.directory, fname)
        if not os.path.exists(robo.directory):
            os.makedirs(robo.directory)
    return full_name


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


class Symoro:
    """Symbol manager, responsible for symbol replacing, file writing."""
    def __init__(self, file_out='disp'):
        """Default values correspond to empty dictionary and screen output.
        """
        self.file_out = file_out
        """Output descriptor. Can be None, 'disp', file
        defines the output destination"""
        self.sydi = {}
        """Dictionary. All the substitutions are saved in it"""
        self.revdi = {}
        """Dictionary. Revers to the self.sydi"""
        self.order_list = []
        """keeps the order of variables to be compute"""

    def simp(self, sym):
        sym = factor(sym)
        new_sym = ONE
        for e in Mul.make_args(sym):
            if e.is_Pow:
                e, p = e.args
            else:
                p = 1
            e = self.C2S2_simp(e)
            e = self.CS12_simp(e)
            new_sym *= e**p
        return new_sym

    def C2S2_simp(self, sym):
        """
        Example
        =======
        >> print C2S2_simp(sympify("-C**2*RL + S*(D - RL*S)"))
        D*S - RL
        """
        if not sym.is_Add:
            repl_dict = {}
            for term in sym.atoms(Add):
                repl_dict[term] = self.C2S2_simp(term)
            sym = sym.xreplace(repl_dict)
            return sym
        names, short_form = trigonometric_info(sym)
        for name in names:
            if short_form:
                C, S = CS_syms(name)
            else:
                C, S = cos(name), sin(name)
            sym = self.try_opt(ONE, None, S**2, C**2, sym)
        return sym

    def CS12_simp(self, sym, silent=False):
        """
        Example
        =======
        >> print Symoro().CS12_simp(sympify("C2*C3 - S2*S3"))
        C23 = C2*C3 - S2*S3
        C23
        >> print Symoro().CS12_simp(sympify("C2*S3*R + S2*C3*R"))
        S23 = C2*S3 + S2*C3
        R*S23
        """
        if not sym.is_Add:
            repl_dict = {}
            for term in sym.atoms(Add):
                repl_dict[term] = self.CS12_simp(term)
            sym = sym.xreplace(repl_dict)
            return sym
        names, short_form = trigonometric_info(sym)
        names = list(names)
        names.sort()
        for n1, n2 in combinations(names, 2):
            if short_form:
                C1, S1 = CS_syms(n1)
                C2, S2 = CS_syms(n2)
                np1, nm1 = get_pos_neg(n1)
                np2, nm2 = get_pos_neg(n2)
                n12 = ang_sum(np1, np2, nm1, nm2)
                nm12 = ang_sum(np1, nm2, nm1, np2)
                C12, S12 = CS_syms(n12)
                C1m2, S1m2 = CS_syms(nm12)
            else:
                C1, S1 = cos(n1), sin(n1)
                C2, S2 = cos(n2), sin(n2)
                C12, S12 = cos(n1+n2), sin(n1+n2)
                C1m2, S1m2 = cos(n1-n2), sin(n1-n2)
            sym = self.try_opt(S12, S1m2, S1*C2, C1*S2, sym, silent)
            sym = self.try_opt(C12, C1m2, C1*C2, -S1*S2, sym, silent)
        return sym

    def try_opt(self, A, Am, B, C, old_sym, silent=False):
        """Replaces B + C by A or B - C by Am.
        Chooses the best option.
        """
        Bcfs = get_max_coef_list(old_sym, B)
        Ccfs = get_max_coef_list(old_sym, C)
        if Bcfs != [] and Ccfs != []:
            Res = old_sym
            Res_tmp = Res
            for coef in Bcfs:
                Res_tmp += A*coef - B*coef - C*coef
                if sym_less(Res_tmp, Res):
                    Res = Res_tmp
            if sym_less(Res, old_sym) and Am is None:
                if not A.is_number and not silent:
                    self.add_to_dict(A, B + C)
                return Res
            elif Am is not None:
                Res2 = old_sym
                Res_tmp = Res2
                for coef in Bcfs:
                    Res_tmp += Am*coef - B*coef + C*coef
                    if sym_less(Res_tmp, Res2):
                        Res2 = Res_tmp
                if sym_less(Res2, Res) and sym_less(Res2, old_sym):
                    if not Am.is_number and not silent:
                        self.add_to_dict(Am, B - C)
                    return Res2
                elif sym_less(Res, old_sym):
                    if not A.is_number and not silent:
                        self.add_to_dict(A, B + C)
                    return Res
        return old_sym

    def add_to_dict(self, new_sym, old_sym):
        """Internal function.
        Extends symbol dictionary by (new_sym, old_sym) pair
        """
        new_sym = sympify(new_sym)
        if new_sym.as_coeff_Mul()[0] == -ONE:
            new_sym = -new_sym
            old_sym = -old_sym
        if new_sym not in self.sydi:
            self.sydi[new_sym] = old_sym
            self.revdi[old_sym] = new_sym
            self.order_list.append(new_sym)
            self.write_equation(new_sym, old_sym)

    def trig_replace(self, M, angle, name):
        """Replaces trigonometric expressions cos(x)
        and sin(x) by CX and SX

        Parameters
        ==========
        M: var or Matrix
            Object of substitution
        angle: var
            symbol that stands for the angle value
        name: int or string
            brief name X for the angle

        Notes
        =====
        The cos(x) and sin(x) will be replaced by CX and SX,
        where X is the name and x is the angle
        """
        if not isinstance(angle, Expr) or angle.is_number:
            return M
        cos_sym, sin_sym = CS_syms(name)
        sym_list = [(cos_sym, cos(angle)), (sin_sym, sin(angle))]
        subs_dict = {}
        for sym, sym_old in sym_list:
            subs_dict[sym_old] = sym
            self.add_to_dict(sym, sym_old)
        for i1 in xrange(M.shape[0]):
            for i2 in xrange(M.shape[1]):
                M[i1, i2] = M[i1, i2].subs(subs_dict)
        return M

    def replace(self, old_sym, name, index='', forced=False):
        """Creates a new symbol for the symbolic expression old_sym.

        Parameters
        ==========
        old_sym: var
            Symbolic expression to be substituted
        name: string or var
            denotion of the expression
        index: int or string, optional
            will be attached to the name. Usualy used for link or joint number.
            Parameter exists for usage convenience
        forced: bool, optional
            If True, the new symbol will be created even if old symbol
            is a simple expression

        Notes
        =====
        Generaly only complex expressions, which contain + - * / ** operations
        will be replaced by a new symbol
        """
        inv_sym = -old_sym
        is_simple = old_sym.is_Atom or inv_sym.is_Atom
        if is_simple and not forced:
            return old_sym
        elif not forced:
            for i in (1, -1):
                if i * old_sym in self.revdi:
                    return i * self.revdi[i * old_sym]
        new_sym = var(str(name) + str(index))
        self.add_to_dict(new_sym, old_sym)
        return new_sym

    def mat_replace(self, M, name, index='',
                    forced=False, skip=0, symmet=False):
        """Replaces each element in M by symbol

        Parameters
        ==========
        M: Matrix
            Object of substitution
        name: string
            denotion of the expression
        index: int or string, optional
            will be attached to the name. Usualy used for link
            or joint number. Parameter exists for usage convenience
        forced: bool, optional
            If True, the new symbol will be created even if old symbol
            is a simple expression
        skip: int, optional
            Number of bottom rows of the matrix, which will be skipped.
            Used in case of Transformation matrix and forced = True.
        symmet: bool, optional
            If true, only for upper triangle part of the matrix
            symbols will be created. The bottom triangle part the
            same symbols will be used


        Returns
        =======
        M: Matrix
            Matrix with all the elements replaced

        Notes
        =====
        -Each element M_ij will be replaced by
            symbol name + i + j + index
        -There are two ways to use this function (examples):
            1)  >>> A = B+C+...
                >>> symo.mat_replace(A, 'A')
                # for the case when expression B+C+... is too big
            2)  >>> A = symo.mat_replace(B+C+..., 'A')
                # for the case when B+C+... is small enough
        """
        for i2 in xrange(M.shape[1]):
            for i1 in xrange(M.shape[0] - skip):
                if symmet and i2 < i1:
                    M[i1, i2] = M[i2, i1]
                    continue
                if M.shape[1] > 1:
                    name_index = name + str(i1 + 1) + str(i2 + 1)
                else:
                    name_index = name + str(i1 + 1)
                M[i1, i2] = self.replace(M[i1, i2], name_index, index, forced)
        return M

    def unfold(self, expr):
        """Unfold the expression using the dictionary.

        Parameters
        ==========
        expr: symbolic expression
            Symbolic expression to be unfolded

        Returns
        =======
        expr: symbolic expression
            Unfolded expression
        """
        while self.sydi.keys() & expr.atoms():
            expr = expr.subs(self.sydi)
        return expr

    def write_param(self, name, header, robo, N):
        """Low-level function for writing the parameters table

        Parameters
        ==========
        name: string
            the name of the table
        header: list
            the table header
        robo: Robot
            Instance of parameter container
        N: list of int
            Indices for which parameter rows will be written
        """
        self.write_line(name)
        self.write_line(l2str(header))
        for j in N:
            params = robo.get_param_vec(header, j)
            self.write_line(l2str(params))
        self.write_line()

    #TDO: rewrite docstring
    def write_params_table(self, robo, title='', geom=True, inert=False,
                           dynam=False, equations=True,
                           inert_name='Dynamic inertia parameters'):
        """Writes the geometric parameters table

        Parameters
        ==========
        robo: Robot
            Instance of the parameter container.
        title: string
            The document title.

        Notes
        =====
        The synamic model generation program can be started with this function
        """
        if title != '':
            self.write_line(title)
            self.write_line()
        if geom:
            self.write_param('Geometric parameters', robo.get_geom_head(),
                             robo, range(1, robo.NF))
        if inert:
            if robo.is_mobile:
                start_frame = 0
            else:
                start_frame = 1

            self.write_param(inert_name, robo.get_dynam_head(),
                             robo, range(start_frame, robo.NL))
        if dynam:
            self.write_param('External forces and joint parameters',
                             robo.get_ext_dynam_head(),
                             robo, range(1, robo.NL))
            self.write_param('Base velicities parameters',
                             robo.get_base_vel_head(),
                             robo, [0, 1, 2])
        if equations:
            self.write_line('Equations:')

    def unknown_sep(self, eq, known):
        """If there is a sum inside trigonometric function and
        the atoms are not the subset of 'known',
        this function will replace the trigonometric symbol bu sum,
        trying to separate known and unknown terms
        """
        while True:
            res = False
            trigs = eq.atoms(sin, cos)
            for trig in trigs:
                args = trig.args[0].atoms()
                if args & known and not args <= known and trig in self.sydi:
                    eq = eq.subs(trig, self.sydi[trig])
                    res = True
            if not res:
                break
        return eq

    def write_equation(self, A, B):
        """Writes the equation A = B into the output

        Parameters
        ==========
        A: expression or var
            left-hand side of the equation.
        B: expression or var
            right-hand side of the equation
        """
        self.write_line(str(A) + ' = ' + str(B))

    def write_line(self, line=''):
        """Writes string data into tha output with new line symbol

        Parameters
        ==========
        line: string, optional
            Data to be written. If empty, it adds an empty line
        """
        if self.file_out == 'disp':
            print line
        elif self.file_out is not None:
            self.file_out.write(str(line) + '\n')

    def file_open(self, robo, ext):
        """
        Initialize file stream

        Parameters
        ==========
        robo: Robot instance
            provides the robot's name
        ext: string
            provides the file name extention
        """
        fname = make_fname(robo, ext)
        self.file_out = open(fname, 'w')

    def file_close(self):
        """
        Initialize file stream

        Parameters
        ==========
        robo: Robot instance
            provides the robot's name
        ext: string
            provides the file name extention
        """
        if self.file_out is not None:
            self.write_line('*=*')
            self.file_out.close()

    def gen_fheader(self, name, *args):
        fun_head = []
        fun_head.append('def %s_func(*args):\n' % name)
        imp_s_1 = 'from numpy import pi, sin, cos, sign\n'
        imp_s_2 = 'from numpy import array, arctan2 as atan2, sqrt\n'
        fun_head.append('    %s' % imp_s_1)
        fun_head.append('    %s' % imp_s_2)
        for i, var_list in enumerate(args):
            v_str_list = self.convert_syms(args[i], True)
            fun_head.append('    %s=args[%s]\n' % (v_str_list, i))
        return fun_head

    def convert_syms(self, syms, rpl_liter=False):
        """Converts 'syms' structure to sintactically correct string

        Parameters
        ==========
        syms: list, Matrix or tuple of them
        rpl_liter: bool
            if true, all literals will be replaced with 'NULx' name.
            It is done to evoid expression like [x, 0] = args[1]
            Because it will cause exception of assigning to literal
        """
        if isinstance(syms, tuple) or isinstance(syms, list):
            syms = [self.convert_syms(item, rpl_liter) for item in syms]
            res = '['
            for i, s in enumerate(syms):
                res += s
                if i < len(syms) - 1:
                    res += ','
            res += ']'
            return res
        elif isinstance(syms, Matrix):
            res = '['
            for i in xrange(syms.shape[0]):
                res += self.convert_syms(list(syms[i, :]), rpl_liter)
                if i < syms.shape[0] - 1:
                    res += ','
            res += ']'
            return res
        else:
            if rpl_liter and sympify(syms).is_number:
                return 'NUL'
            else:
                return str(syms)

    def extract_syms(self, syms):
        """ returns set of all symbols from list or matrix
        or tuple of them
        """
        if isinstance(syms, tuple) or isinstance(syms, list):
            atoms = (self.extract_syms(item) for item in syms)
            return reduce(set.__or__, atoms, set())
        elif isinstance(syms, Matrix):
            return self.extract_syms(list(syms))
        elif isinstance(syms, Expr):
            return syms.atoms(Symbol)
        else:
            return set()

    def sift_syms(self, rq_syms, wr_syms):
        """Returns ordered list of variables to be compute
        """
        order_list = []   # vars that are defined in sydi
        for s in reversed(self.order_list):
            if s in rq_syms and not s in wr_syms:
                order_list.insert(0, s)
                s_val = self.sydi[s]
                if isinstance(s_val, Expr):
                    atoms = s_val.atoms(Symbol)
                    rq_syms |= {s for s in atoms if not s.is_number}
        rq_vals = [s for s in rq_syms if not (s in self.sydi or s in wr_syms)]
            # required vars that are not defined in sydi
            # will be set to '1.'
        return rq_vals + order_list

    def gen_fbody(self, name, to_return, wr_syms, multival):
        """Generates list of string statements of the function that
        computes symbolf from to_return.  wr_syms are considered to
        be known
        """
        # final symbols to be compute
        syms = self.extract_syms(to_return)
        # defines order of computation
        order_list = self.sift_syms(syms, wr_syms)
        # list of instructions in final function
        fun_body = []
        # will be switched to true when branching detected
        space = '    '
        folded = 1    # indentation = 1 + number of 'for' statements

        for s in order_list:
            if s not in self.sydi:
                item = '%s%s=1.\n' % (space * folded, s)
            elif isinstance(self.sydi[s], tuple):
                multival = True
                item = '%sfor %s in %s:\n' % (space * folded, s, self.sydi[s])
                folded += 1
            else:
                item = '%s%s=%s\n' % (space * folded, s, self.sydi[s])
            fun_body.append(item)
        ret_expr = self.convert_syms(to_return)
        if multival:
            fun_body.insert(0, '    %s_result=[]\n' % (name))
            item = '%s%s_result.append(%s)\n' % (space*folded, name, ret_expr)
        else:
            item = '    %s_result=%s\n' % (name, ret_expr)
        fun_body.append(item)
        fun_body.append('    return %s_result\n' % (name))
        return fun_body

    def gen_func(self, name, to_return, args, multival=False):
        """ Returns function that computes what is in to_return
        using *args as arguments

         Parameters
        ==========
        name: string
            Future function's name, must be different for
            different fucntions
        to_return: list, Matrix or tuple of them
            Determins the shape of the output and symbols inside it
        *args: any number of lists, Matrices or tuples of them
            Determins the shape of the input and symbols
            names to assigned

        Notes
        =====
        -All unassigned used symbols will be set to '1.0'.
        -This function must be called only after the model that
            computes symbols in to_return have been generated.
        """
        fun_head = self.gen_fheader(name, args)
        wr_syms = self.extract_syms(args)   # set of defined symbols
        fun_body = self.gen_fbody(name, to_return, wr_syms, multival)
        fun_string = "".join(fun_head + fun_body)
        exec fun_string
        print fun_string
#  TODO:       print is for debug pupuses, to be removed
        return eval('%s_func' % name)
