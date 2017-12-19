# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package provides description
of the robot parametrizaion container and symbol replacer class.
"""


import os
import copy

from sympy import sin, sign
from sympy import Matrix, Expr
from sympy import zeros, var, sympify, eye

from . import baseparams
from . import dyniden
from . import inertia
from . import nealgos
from symoroutils import filemgr
from symoroutils import symbolmgr
from symoroutils.tools import FAIL, OK
from symoroutils.tools import TREE, INT_KEYS


class Robot(object):
    """Container of the robot parametric description.

    Responsible for low-level geometric transformation
    and direct geometric model generation.
    Also provides different representations of parameters.
    """
    def __init__(
        self, name, NL=0, NJ=0, NF=0, is_floating=False,
        structure=TREE, is_mobile=False, directory=None,
        par_file_path=None
    ):
        # Name of the robot: string.
        self.name = name
        # Directory name.
        self.directory = self.set_directory(directory)
        # PAR file path.
        self.par_file_path = self.set_par_file_path(par_file_path)
        # Whether the base frame is floating: bool.
        self.is_floating = is_floating
        # Whether the robot is a mobile robot: bool.
        self.is_mobile = is_mobile
        # Number of links: int.
        self.nl = NL
        # Number of joints: int.
        self.nj = NJ
        # Number of frames: int.
        self.nf = NF
        # Type of robot structure.
        self.structure = structure
        # Joint type: list of int.
        self.sigma = [0] * (NF + 1)
        # Index of antecedent joint: list of int.
        self.ant = list(range(-1, self.NF - 1))
        # Actuated, if 1, then the joint is actuated.
        self.mu = [1] * (NF + 1)
        # Geometrical parameter: list of var.
        self.theta = [0] + [var('th{}'.format(i + 1)) for i in range(NF)]
        # Geometrical parameter: list of var.
        self.r = [0] * (NF + 1)
        # Geometrical parameter: list of var.
        self.alpha = [0] * (NF + 1)
        # Geometrical parameter: list of var.
        self.d = [0] * (NF + 1)
        # Geometrical parameter: list of var.
        self.gamma = [0] * (NF + 1)
        # Geometrical parameter: list of var.
        self.b = [0] * (NF + 1)
        # Transformation from reference frame to zero frame.
        self.Z = eye(4)
        num = list(range(self.NL))
        numj = list(range(self.NJ))
        # Base angular velocity: 3x1 matrix.
        self.w0 = zeros(3, 1)
        # Base angular acceleration: 3x1 matrix.
        self.wdot0 = zeros(3, 1)
        # Base linear velocity: 3x1 matrix.
        self.v0 = zeros(3, 1)
        # Base linear acceleration: 3x1 matrix.
        self.vdot0 = zeros(3, 1)
        # Joint speed: list of var.
        self.qdot = [var('QP{0}'.format(i)) for i in numj]
        # Joint acceleration: list of var.
        self.qddot = [var('QDP{0}'.format(i)) for i in numj]
        # External moment of link: list of 3x1 matrix.
        self.Nex = [zeros(3, 1) for i in num]
        self.Nex[-1] = Matrix(var('CX{0}, CY{0}, CZ{0}'.format(self.NL - 1)))
        # External force of link: list of 3x1 matrix.
        self.Fex = [zeros(3, 1) for i in num]
        self.Fex[-1] = Matrix(var('FX{0}, FY{0}, FZ{0}'.format(self.NL - 1)))
        # Dry friction coefficient: list of var.
        self.FS = [var('FS{}'.format(i)) for i in num]
        # Joint actuator inertia: list of var.
        self.IA = [var('IA{}'.format(i)) for i in num]
        # Viscous friction coefficient: list of var.
        self.FV = [var('FV{}'.format(i)) for i in num]
        # First momentum of link: list of 3x1 matrix.
        self.MS = [Matrix(var('MX{0}, MY{0}, MZ{0}'.format(i))) for i in num]
        # Mass of link: list of var.
        self.M = [var('M{}'.format(i)) for i in num]
        # Joint torques: list of var.
        self.GAM = [var('GAM{}'.format(i)) for i in numj]
        # Inertia tensor of link: list of 3x3 matrix.
        J_str = 'XX{0},XY{0},XZ{0},XY{0},YY{0},YZ{0},XZ{0},YZ{0},ZZ{0}'
        self.J = [Matrix(3, 3, var(J_str.format(i))) for i in num]
        # Gravity vector: 3x1 matrix.
        self.G = Matrix([0, 0, var('GZ')])
        # eta - rigid or flexible.
        self.eta = [0] * self.NJ
        # k - joint stiffness.
        self.k = [0] * self.NJ

    def set_par_file_path(self, path=None):
        if path is None or not os.path.isabs(path):
            file_path = filemgr.get_file_path(self)
        else:
            file_path = path
        self.par_file_path = file_path
        return file_path

    def set_directory(self, path=None):
        if path is None or not os.path.isdir(path):
            directory = filemgr.get_folder_path(self.name)
        else:
            directory = path
        self.directory = directory
        return directory

    @property
    def par_file_name(self):
        """Return the PAR file name"""
        head, tail = os.path.split(self.par_file_path)
        return tail.strip()

    def set_defaults(self, joint=False, geom=False, base=False):
        # Joint parameters.
        if joint:
            self._set_joint_defaults()
        # Geometric parameters.
        if geom:
            # do nothing
            #self._set_geom_defaults()
            pass
        # Base parameters.
        if base:
            self._set_base_defaults()

    def put_val(self, j, name, val):
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
        elif name == 'eta':
            self.eta[j] = int(val)
        elif name == 'k':
            self.k[j] = val
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
        elif name == 'eta':
            return self.eta[j]
        elif name == 'k':
            return self.k[j]
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

    def compute_idym(self):
        """
        Compute the Inverse Dynamic Model of the robot using the
        recursive Newton-Euler algorithm. Also choose the Newton-Euler
        algorithm based on the robot type.
        """
        symo = symbolmgr.SymbolManager()
        symo.file_open(self, 'idm')
        title = "Inverse Dynamic Model using Newton-Euler Algorithm\n"
        if 1 in self.eta:
            # with flexible joints
            title = title + "Robot with flexible joints\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.flexible_inverse_dynmodel(self, symo)
        elif self.is_floating:
            # with rigid joints and floating base
            title = title + "Robot with rigid joints and floating base\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.composite_inverse_dynmodel(self, symo)
        elif self.is_mobile:
            # mobile robot with rigid joints - known base acceleration
            title = title + "Robot with mobile base (Vdot0 is known)\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.mobile_inverse_dynmodel(self, symo)
        else:
            # with rigid joints and fixed base
            title = title + "Robot with rigid joints and fixed base\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.fixed_inverse_dynmodel(self, symo)
        symo.file_close()
        return symo

    def compute_inertiamatrix(self):
        """
        Compute the Inertia Matrix of the robot using the Composite link
        algorithm.
        """
        symo = symbolmgr.SymbolManager()
        symo.file_open(self, 'inm')
        title = "Inertia matrix using Composite links algorithm\n"
        if self.is_floating or self.is_mobile:
            # with floating or mobile base
            title = title + "Robot with floating/mobile base\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            inertia.floating_inertia_matrix(self, symo)
        else:
            # with fixed base
            title = title + "Robot with fixed base\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            inertia.fixed_inertia_matrix(self, symo)
        symo.file_close()
        return symo

    def compute_ddym(self):
        """
        Compute the Direct Dynamic Model of the robot using the
        recursive Newton-Euler algorithm.
        """
        symo = symbolmgr.SymbolManager()
        symo.file_open(self, 'ddm')
        title = "Direct Dynamic Model using Newton-Euler Algorithm\n"
        if self.is_floating:
            # with floating base
            title = title + "Robot with floating base\n"
        else:
            # with fixed base
            title = title + "Robot with fixed base\n"
        symo.write_params_table(self, title, inert=True, dynam=True)
        nealgos.direct_dynmodel(self, symo)
        symo.file_close()
        return symo

    def compute_pseudotorques(self):
        """
        Compute Coriolis, Centrifugal, Gravity, Friction and external
        torques using Newton-Euler algortihm.
        """
        pseudo_robo = copy.deepcopy(self)
        pseudo_robo.qddot = zeros(pseudo_robo.NL, 1)
        symo = symbolmgr.SymbolManager()
        symo.file_open(self, 'ccg')
        title = "Pseudo forces using Newton-Euler Algorithm\n"
        if 1 in pseudo_robo.eta:
            # with flexible joints
            title = title + "Robot with flexible joints\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.flexible_inverse_dynmodel(pseudo_robo, symo)
        elif pseudo_robo.is_floating:
            # with rigid joints and floating base
            title = title + "Robot with rigid joints and floating base\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.composite_inverse_dynmodel(pseudo_robo, symo)
        elif pseudo_robo.is_mobile:
            # mobile robot with rigid joints - known base acceleration
            title = title + "Robot with mobile base (Vdot0 is known)\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.mobile_inverse_dynmodel(pseudo_robo, symo)
        else:
            # with rigid joints and fixed base
            title = title + "Robot with rigid joints and fixed base\n"
            symo.write_params_table(self, title, inert=True, dynam=True)
            nealgos.fixed_inverse_dynmodel(pseudo_robo, symo)
        symo.file_close()
        return symo

    def compute_baseparams(self):
        """
        Compute the Base Inertial Parameters of the robot.
        """
        base_robo = copy.deepcopy(self)
        symo = symbolmgr.SymbolManager()
        symo.file_open(base_robo, 'regp')
        title = "Base Inertial Parameters equations"
        symo.write_params_table(
            base_robo, title, inert=True, dynam=True
        )
        # compute base inertial params
        baseparams.base_inertial_parameters(base_robo, symo)
        symo.write_line()
        title = "Grouped inertial parameters"
        symo.write_params_table(
            base_robo, title, inert=True, equations=False
        )
        symo.file_close()
        # set new name for robot with base params
        base_robo.name = base_robo.name + "_base"
        file_path = filemgr.get_file_path(base_robo)
        base_robo.set_par_file_path(file_path)
        return symo, base_robo

    def compute_dynidenmodel(self):
        """
        Compute the Dynamic Identification model of the robot.
        """
        symo = symbolmgr.SymbolManager()
        symo.file_open(self, 'dim')
        title = "Dynamic Identification Model (Newton-Euler method)"
        symo.write_params_table(self, title, inert=True, dynam=True)
        dyniden.dynamic_identification_model(self, symo)
        symo.file_close()
        return symo

    @property
    def q_vec(self):
        """Generates vector of joint variables
        """
        qs = []
        for i in range(1, self.NJ):
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
        for i in range(1, self.NJ):
            if self.mu[i] == 0:
                q.append(self.get_q(i))
        return q

    @property
    def q_active(self):
        """Generates vector of active joint variables (including cut!)
        """
        q = list()
        for i in range(1, self.NJ):
            if self.mu[i] == 1:
                q.append(self.get_q(i))
        return q

    @property
    def indx_passive(self):
        """Generates vector of passive joint indices
        """
        return [i for i in range(1, self.NL) if self.mu[i] == 0]

    @property
    def indx_active(self):
        """Generates vector of active joint indices
        """
        return [i for i in range(1, self.NL) if self.mu[i] == 1]

    @property
    def indx_cut(self):
        """Generates vector of cut joint indices
        """
        return list(range(self.NL, self.NJ))

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
        return self.v0

    @property
    def VP0(self):
        return self.vdot0

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
        return [(i, i + B) for i in range(self.NL, self.NJ)]

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
        j: int, Joint index.

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
        if j not in range(self.NF):
            return angs
        if type(self.theta[j]) != int and not self.theta[j].is_number:
            angs.append((self.theta[j], j))
        if type(self.alpha[j]) != int and not self.alpha[j].is_number:
            angs.append((self.alpha[j], 'A{}'.format(j)))
        if type(self.gamma[j]) != int and not self.gamma[j].is_number:
            angs.append((self.gamma[j], 'G{}'.format(j)))
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
        """Returns 10-vector of inertia parameters of link j.

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
                'FS', 'FV', 'QP', 'QDP', 'GAM', 'eta', 'k']

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

    def _set_joint_defaults(self):
        """
        Set default values for joint parameters for those exceptional
        from the ones set in the ctor.
        """
        for j in range(1, self.NJ):
            try:
                if self.sigma[j] == 2:
                    self.qdot[j] = 0
                    self.qddot[j] = 0
                    self.GAM[j] = 0
                else:
                    self.qdot[j] = var('QP{0}'.format(j))
                    self.qddot[j] = var('QDP{0}'.format(j))
                    self.GAM[j] = var('GAM{0}'.format(j))
            except IndexError:
                # just ignore exception
                pass
            if self.eta[j] == 1:
                self.k[j] = var('k{0}'.format(j))
            else:
                self.k[j] = 0

    def _set_geom_defaults(self):
        """
        Set default values for geometric parameters for those
        exceptional from the ones set in the ctor.
        """
        for j in range(1, self.NF):
            if self.sigma[j] == 0:
                self.theta[j] = var('th{0}'.format(j))
            elif self.sigma[j] == 1:
                self.r[j] = var('r{0}'.format(j))
            elif self.sigma[j] == 2:
                self.mu[j] = 0

    def _set_base_defaults(self):
        """
        Set default values for base parameters for those exceptional
        from the ones set in the ctor.
        """
        if self.is_floating or self.is_mobile:
            self.G = Matrix([var('GX'), var('GY'), var('GZ')])
            self.v0 = Matrix([var('VXb'), var('VYb'), var('VZb')])
            self.w0 = Matrix([var('WXb'), var('WYb'), var('WZb')])
            self.vdot0 = Matrix([var('VPXb'), var('VPYb'), var('VPZb')])
            self.wdot0 = Matrix([var('WPXb'), var('WPYb'), var('WPZb')])
            # Z matrix
            for i in range(0, 3):
                for j in range(0, 3):
                    self.Z[i, j] = var('Zr{0}{1}'.format(i + 1, j + 1))
            for j in range(0, 3):
                self.Z[j, 3] = var('Zt{0}'.format(j + 1))
