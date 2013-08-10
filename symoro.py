"""
This module of SYMORO package provides description of the robot parametrizaion
container and symbol replacer class.

The core symbolic library is sympy.

ECN - ARIA1 2013
"""
from sympy import Matrix,zeros,var,sympify,pi,sign, Integer
from copy import copy
import re

ZERO = Integer(0)

class Robot:
    """Container of the robot parametric description.
    Responsible for low-level geometric transformation
    and direct geometric model generation.
    Also provides different representations of parameters."""

    #member variables:
    name = 'Empty'
    """  name of the robot : string"""
    NL = 0
    """  number of links : int"""
    NJ = 0
    """  number of joints : int"""
    NF = 0
    """  number of frames : int"""
    sigma = None
    """  joint type : list of int"""
    ant = None
    """  index of antecedent joint : list of int"""
    num = None
    """order numbers of joints (for display purposes) : list of int.
                The last number corresponds to the base frame -1 """
    mu = None
    """motorization, if 1, then the joint im motorized"""
    theta = None
    """  geometrical parameter : list of var"""
    r = None
    """  geometrical parameter  : list of var"""
    alpha = None
    """  geometrical parameter  : list of var"""
    d = None
    """  geometrical parameter : list of var"""
    gamma = None
    """  geometrical parameter : list of var"""
    b = None
    """  geometrical parameter : list of var"""
    J = None
    """  inertia tensor of link : list of 3x3 matrix"""
    MS = None
    """  first momentum of link : list of 3x1 matrix"""
    M = None
    """  mass of link : list of var"""
    G = (0,0,0)
    """  gravity vector : 3x1 matrix"""
    GAM = None
    """  joint torques : list of var"""
    w0 = (0,0,0)
    """  base angular velocity : 3x1 matrix"""
    wdot0 = (0,0,0)
    """  base angular acceleration : 3x1 matrix"""
    v0 = (0,0,0)
    """  base linear velocity : 3x1 matrix"""
    vdot0 = (0,0,0)
    """  base linear acceleration : 3x1 matrix"""
    qdot = None
    """  joint speed : list of var"""
    qddot = None
    """  joint acceleration : list of var"""
    Nex = None
    """  external moment of link : list of 3x1 matrix"""
    Fex = None
    """  external force of link : list of 3x1 matrix"""
    FS = None
    """  dry friction coefficient : list of ver"""
    FV = None
    """  fluid friction coefficient : list of var"""
    IA = None
    """  joint actuator inertia : list of var"""

    #member methods:
    def gen_q_vec(self):
        """Generates vector of joint variables
        """
        self.q = []
        for i in range(self.NL):
            if self.sigma[i] == 0:
                self.q.append(self.theta[i])
            elif self.sigma[i] == 1:
                self.q.append(self.r[i])
            else:
                self.q.append(0)

    def fric_v(self,j):
        """Fluid friction torque

        Parameters
        ==========
        j : int
            Joint index.

        Returns
        =======
        fric_v : sympy expression
            Expression for fluid friction torque of joint j
        """
        return self.FV[j]*self.qdot[j]

    def fric_s(self,j):
        """Dry friction torque

        Parameters
        ==========
        j : int
            Joint index.

        Returns
        =======
        fric_s : sympy expression
            Expression for dry friction torque of joint j
        """
        return self.FS[j]*sign(self.qdot[j])

    def tau_ia(self,j):
        """Actuator inertia torque

        Parameters
        ==========
        j : int
            Joint index.

        Returns
        =======
        fric_v : sympy expression
            Expression for actuator inertia torque of joint j
        """
        return self.IA[j]*self.qddot[j]

    def get_angles(self,j):
        """List of non-constant angles of frame j

        Parameters
        ==========
        j : int
            Frame index.

        Returns
        =======
        get_angles : list of touples (var,name)
            Returns list of touples, where:
            var - the angle symbol,
            name - brief name for cos and sin abbreviation
        """
        angs = []
        if j not in range(self.NF):
            return angs
        index = str(self.num[j])
        if not sympify(self.theta[j]).is_constant():
            angs.append((self.theta[j],index))
        if not sympify(self.alpha[j]).is_constant():
            angs.append((self.alpha[j],'A' + index))
        if not sympify(self.gamma[j]).is_constant():
            angs.append((self.gamma[j],'G' + index))
        return angs

    def chain(self,j,k = - 1):
        """Chain of antecedent frames between j-th and k-th frames

        Parameters
        ==========
        j : int
            Start frame index.
        k : int
            Final frame index.

        Returns
        =======
        u : list of ints
            List of antecedent frames. j is the first index in the list.
            k is not included
        """
        u = []
        while j != k and j != -1:
            u.append(j)
            j = self.ant[j]
        return u

    def common_root(self,i,j):
        """Common root j-th and i-th frames

        Parameters
        ==========
        j : int
            Frame index.
        i : int
            Frame index.

        Returns
        =======
        common_root : int
            The highest index of the common frame in chains for i and j.
            If they don't have common root, -1
        """
        u = self.chain(i)
        while j != - 1:
            if j in u:
                return j
            j = self.ant[j]
        return  - 1

    def put_dynam_param(self,K,j):
        """Write the inertia parameters of link j from 10-vector K.

        Parameters
        ==========
        K : Matrix 10x1
            Vector of inertia parameters
        j : int
            Link index.
        """
        self.J[j] = Matrix([[K[0],K[1],K[2]],
                    [K[1],K[3],K[4]],
                    [K[2],K[4],K[5]]])
        self.MS[j] = Matrix(3,1,K[6:9])
        self.M[j] = K[9]

    def get_ext_dynam_head(self):
        """Returns header for external forces and torques, friction parameters
        and joint speeds, accelerations. Used for output generation.

        Returns
        =======
        get_ext_dynam_head : list of strings
        """
        return ['j','FX','FY','FZ','CX','CY','CZ','FS','FV','QP','QDP','GAM']

    def get_inert_head(self):
        """Returns header for inertia parameters. Used for output generation.

        Returns
        =======
        get_inert_head : list of strings
        """
        return ['j','XX','XY','XZ','YY','YZ','ZZ','MX','MY','MZ','M','IA']

    def get_geom_head(self):
        """Returns header for geometric parameters. Used for output generation.

        Returns
        =======
        get_geom_head : list of strings
        """
        return ['j','ant','sigma','gamma','b','alpha','d','theta','r']

    def get_base_vel_head(self):
        """Returns header for base velocities and gravity vector.
        Used for output generation.

        Returns
        =======
        get_base_vel_head : list of strings
        """
        return ['j','W0','WP0','V0','VP0','G']

    def get_geom_param(self,j):
        """Returns vector of geometric parameters of frame j.
        Used for output generation.

        Parameters
        ==========
        j : int
            Frame index.

        Returns
        =======
        params : list
        """
        params = [self.num[j],self.num[self.ant[j]],self.sigma[j],
                  self.gamma[j],self.b[j],self.alpha[j],self.d[j],
                  self.theta[j],self.r[j]]
        return params

    def get_ext_dynam_param(self,j):
        """Returns vector of external forces and torques, friction parameters
        and joint speeds, accelerations of link j.
        Used for output generation.

        Parameters
        ==========
        j : int
            Link index.

        Returns
        =======
        params : list
        """
        params = [self.num[j],self.Fex[j][0],self.Fex[j][1],self.Fex[j][2],
                  self.Nex[j][0],self.Nex[j][1],self.Nex[j][2],
                  self.FS[j],self.FV[j],self.qdot[j],
                  self.qddot[j],self.GAM[j]]
        return params

    def get_base_vel(self,j):
        """Returns vector of j-th components of base
        velocities and gravity vector.
        Used for output generation.

        Parameters
        ==========
        j : int
            Link index.

        Returns
        =======
        params : list
        """
        params = [j + 1,self.w0[j],self.wdot0[j],self.v0[j],
                  self.vdot0[j],self.G[j]]
        return params

    def get_inert_param(self,j):
        """Returns vector of inertia paremeters of link j.
        Used for output generation.

        Parameters
        ==========
        j : int
            Link index.

        Returns
        =======
        params : list
        """
        params = [self.num[j],self.J[j][0],self.J[j][1],self.J[j][2],
                  self.J[j][4],self.J[j][5],self.J[j][8],self.MS[j][0],
                  self.MS[j][1],self.MS[j][2],self.M[j],self.IA[j]]
        return params

    def get_dynam_param(self,j):
        """Returns 10-vector of inertia paremeters of link j.

        Parameters
        ==========
        j : int
            Link index.

        Returns
        =======
        get_dynam_param : Matrix 10x1
        """
        K = [self.J[j][0],self.J[j][1],self.J[j][2],self.J[j][4],
                    self.J[j][5],self.J[j][8],self.MS[j][0],self.MS[j][1],
                    self.MS[j][2],self.M[j]]
        return Matrix(K)

    @classmethod
    def CartPole(self):
        """Generates Robot instance of classical CartPole dynamic system."""
        robo = Robot()
        robo.name = 'CartPole'
        robo.ant = ( - 1,0)
        robo.sigma = (1,0)
        robo.alpha = (pi/2,pi/2)
        robo.d = (0,0)
        robo.theta = (pi/2,var('Th2'))
        robo.r = (var('R1'),0)
        robo.b = (0,0)
        robo.gamma = (0,0)
        robo.num = range(1,3)
        robo.NJ = 2
        robo.NL = 2
        robo.NF = 2
        robo.gen_q_vec()
        robo.Nex = [zeros(3,1) for i in robo.num]
        robo.Fex = [zeros(3,1) for i in robo.num]
        robo.FS = [0 for i in robo.num]
        robo.IA = [0 for i in robo.num]
        robo.FV = [var('FV{0}'.format(i)) for i in robo.num]
        robo.MS = [zeros(3,1) for i in robo.num]
        robo.MS[1][0] = var('MX2')
        robo.M = [var('M{0}'.format(i)) for i in robo.num]
        robo.GAM = [var('GAM{0}'.format(i)) for i in robo.num]
        robo.J = [zeros(3) for i in robo.num]
        robo.J[1][2,2] = var('ZZ2')
        robo.G = Matrix([0,0, - var('G3')])
        robo.w0 = zeros(3,1)
        robo.wdot0 = zeros(3,1)
        robo.v0 = zeros(3,1)
        robo.vdot0 = zeros(3,1)
        robo.q = var('R1,Th2')
        robo.qdot = var('R1d,Th2d')
        robo.qddot = var('R1dd,Th2dd')
        robo.num.append(0)
        return robo

    @classmethod
    def SR400(self):
        """Generates Robot instance of RX90"""
        robo = Robot()
        #table of geometric parameters RX90
        robo.name = 'SR400'
        robo.NJ = 9
        robo.NL = 8
        robo.NF = 10
        robo.num = range(1,robo.NF + 1)
        numL = range(1,robo.NL + 1)
        robo.ant = (-1,0,1,2,3,4,0,6,7,2)
        robo.sigma = (0,0,0,0,0,0,0,0,0,2)
        robo.mu = (1,1,0,1,1,1,1,0,0,0)
        robo.alpha = (0,-pi/2,0, - pi/2,pi/2, - pi/2, - pi/2,0,0,0)
        d_var = var('D:9')
        robo.d = (0,d_var[2],d_var[2],d_var[2],0,0,d_var[2],d_var[8],d_var[3],d_var[8])
        robo.theta = list(var('th1:10'))+[0]
        robo.r = (0,0,0,var('RL4'),0,0,0,0,0,0)
        robo.b = (0,0,0,0,0,0,0,0,0,0)
        robo.gamma = (0,0,0,0,0,0,0,0,0,pi/2)
        robo.w0 = zeros(3,1)
        robo.wdot0 = zeros(3,1)
        robo.v0 = zeros(3,1)
        robo.vdot0 = zeros(3,1)
        robo.gen_q_vec()
        robo.qdot = [var('QP{0}'.format(i)) for i in numL]
        robo.qddot = [var('QDP{0}'.format(i)) for i in numL]
        robo.Nex= [zeros(3,1) for i in numL]
        robo.Nex[ - 1] = Matrix(var('CX{0},CY{0},CZ{0}'.format(robo.NL)))
        robo.Fex = [zeros(3,1) for i in numL]
        robo.Fex[ - 1] = Matrix(var('FX{0},FY{0},FZ{0}'.format(robo.NL)))
        robo.FS = [var('FS{0}'.format(i)) for i in numL]
        robo.IA = [var('IA{0}'.format(i)) for i in numL]
        robo.FV = [var('FV{0}'.format(i)) for i in numL]
        robo.MS = [Matrix(var('MX{0},MY{0},MZ{0}'.format(i))) for i in numL]
        robo.M = [var('M{0}'.format(i)) for i in numL]
        robo.GAM = [var('GAM{0}'.format(i)) for i in numL]
        robo.J = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                            'XY{0},YY{0},YZ{0},'
                            'XZ{0},YZ{0},ZZ{0}').format(i))) for i in numL]
        robo.G = Matrix([0,0,var('G3')])
        robo.num.append(0)
        return robo

    @classmethod
    def RX90(self):
        """Generates Robot instance of RX90"""
        robo = Robot()
        #table of geometric parameters RX90
        robo.name = 'RX90'
        robo.NJ = 6
        robo.NL = 6
        robo.NF = 6
        robo.num = range(1,robo.NJ + 1)
        robo.ant = range( - 1,robo.NJ - 1)
        robo.sigma = (0,0,0,0,0,0)
        robo.alpha = (0,pi/2,0, - pi/2,pi/2, - pi/2)
        robo.d = (0,0,var('D3'),0,0,0)
        robo.theta = list(var('th1:7'))
        robo.r = (0,0,0,var('RL4'),0,0)
        robo.b = (0,0,0,0,0,0)
        robo.gamma = (0,0,0,0,0,0)
        robo.w0 = zeros(3,1)
        robo.wdot0 = zeros(3,1)
        robo.v0 = zeros(3,1)
        robo.vdot0 = zeros(3,1)
        robo.gen_q_vec()
        robo.qdot = [var('QP{0}'.format(i)) for i in robo.num]
        robo.qddot = [var('QDP{0}'.format(i)) for i in robo.num]
        robo.Nex= [zeros(3,1) for i in robo.num]
        robo.Nex[ - 1] = Matrix(var('CX{0},CY{0},CZ{0}'.format(robo.num[ - 1])))
        robo.Fex = [zeros(3,1) for i in robo.num]
        robo.Fex[ - 1] = Matrix(var('FX{0},FY{0},FZ{0}'.format(robo.num[ - 1])))
        robo.FS = [var('FS{0}'.format(i)) for i in robo.num]
        robo.IA = [var('IA{0}'.format(i)) for i in robo.num]
        robo.FV = [var('FV{0}'.format(i)) for i in robo.num]

        robo.MS = [Matrix(var('MX{0},MY{0},MZ{0}'.format(i))) for i in robo.num]
        robo.M = [var('M{0}'.format(i)) for i in robo.num]
        robo.GAM = [var('GAM{0}'.format(i)) for i in robo.num]
        robo.J = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                            'XY{0},YY{0},YZ{0},'
                            'XZ{0},YZ{0},ZZ{0}').format(i))) for i in robo.num]
        robo.G = Matrix([0,0,var('G3')])
        robo.num.append(0)
        return robo

class Init:
    @classmethod
    def init_Jplus(self,robo):
        """Copies the inertia parameters. Used for composed link inertia computation

        Returns
        =======
        Jplus: list of Matrices 3x3
        MSplus: list of Matrices 3x1
        Mplus : list of var
        """
        Jplus = copy(robo.J)
        Jplus.append(zeros(3,3))
        MSplus = copy(robo.MS)
        MSplus.append(zeros(3,1))
        Mplus = copy(robo.M)
        Mplus.append(0)
        return Jplus,MSplus,Mplus

    @classmethod
    def init_mat(self,robo,N = 3):
        """Generates a list of Matrices.Size of the list is number of links.

        Parameters
        ==========
        robo : Robot
            Instance of robot description container
        N : int, optional
            size of the matries, default is 3

        Returns
        =======
        list of Matrices NxN
        """
        return [zeros(N,N) for i in range(robo.NL)]

    @classmethod
    def init_vec(self,robo,N = 3,ext=0):
        """Generates a list of vectors. Size of the list is number of links.

        Parameters
        ==========
        robo : Robot
            Instance of robot description container
        N : int, optional
            size of the vectors, default is 3
        ext : int, optional
            additional vector instances over number of links

        Returns
        =======
        list of Matrices Nx1
        """
        return [zeros(N,1) for i in range(robo.NL+ext)]

    @classmethod
    def init_scalar(self,robo):
        """Generates a list of vars. Size of the list is number of links."""
        return [0 for i in range(robo.NL)]

    @classmethod
    def init_w(self,robo):
        """Generates a list of vectors for angular velocities.
        Size of the list is number of links + 1.
        The last vector is the base angular velocity"""
        w = self.init_vec(robo)
        w.append(robo.w0)
        return w

    @classmethod
    def init_wv_dot(self,robo):
        """Generates lists of vectors for angular and linear accelerations.
        Size of the list is number of links + 1.
        The last vector is the base angular velocity

        Returns
        =======
        vdot : list of Matrices 3x1
        wdot : list of Matrices 3x1
        """
        wdot = self.init_vec(robo)
        wdot.append(robo.wdot0)
        vdot = self.init_vec(robo)
        vdot.append(robo.vdot0 - robo.G)
        return wdot,vdot

    @classmethod
    def init_U(self,robo):
        """Generates a list of auxiliary U matrices"""
        U = Init.init_mat(robo)
        #the value for the -1th base frame
        U.append(hat(robo.w0)**2 + hat(robo.wdot0))
        return U

    @classmethod
    def product_combinations(self,v):
        """Generates 6-vector of different v elements' product combinations

        Parameters
        ==========
        v : Matrix 3x1
            vector

        Returns
        =======
        product_combinations : Matrix 6x1
        """
        return Matrix([v[0]*v[0],v[0]*v[1],v[0]*v[2],
                     v[1]*v[1],v[1]*v[2],v[2]*v[2]])

def hat(v):
    """Generates vectorial preproduct matrix

    Parameters
    ==========
    v : Matrix 3x1
        vector

    Returns
    =======
    hat : Matrix 3x3
    """
    return Matrix([[0, - v[2],v[1]],
                   [v[2],0, - v[0]],
                   [ - v[1],v[0],0]])


def l2str(l,spacing = 7):
    """Converts a list into string, that will be written into the text table.

    Parameters
    ==========
    l : list
        List to be converted
    spacing : int, optional
        Defines the size of one cell of the table

    Returns
    =======
    s : string
        String representation

    Notes
    =====
    l2str([1,2,3]) will be converted into '1      2      3      '
    """
    s = ''
    for i in l:
        s += str(i) + ' '*(spacing-len(str(i)))
    return s

def get_trig_couple_names(sym):
    names_s = find_trig_names(sym, r'S', 1)
    names_c = find_trig_names(sym, r'C', 1)
    return  names_c, names_s


def find_trig_names(sym, pref = r'', pref_len = 0, post = r'', post_len = 0):
    search_res = re.findall(pref + r'[AGm0-9]*'+ post,str(sym))
    if post_len == 0:
        return set([s[pref_len:] for s in search_res])
    else:
        return set([s[pref_len:-post_len] for s in search_res])

def get_trig_pow_names(sym, min_pow = 2):
    post = r'\*\*[{0}-9]'.format(min_pow)
    names_s = find_trig_names(sym, r'S', 1, post, 3 )
    names_c = find_trig_names(sym, r'C', 1, post, 3 )
    return names_c & names_s

def get_max_coef(sym, X):
    EC = ZERO
    for EC_cand in sym.as_ordered_terms():
        if (EC_cand/X).as_numer_denom()[1] == 1:
            EC += EC_cand/X
    return EC

def ang_sum(np1,np2,nm1,nm2):
    np2,nm1 = reduce_str(np2,nm1)
    np1,nm2 = reduce_str(np1,nm2)
    if len(nm1) + len(nm2) == 0: return np1 + np2
    else: return np1 + np2 + 'm' + nm1 + nm2

def get_pos_neg(s):
    if s.find('m') != -1:
        return s.split('m')[0],s.split('m')[1]
    else:
        return s,''

def reduce_str(s1,s2):
    for j,char in enumerate(s1):
        i = s2.find(char)
        if s2.find(char) != -1:
            s2 = s2[:i] + s2[i+1:]
            s1 = s1[:j] + s1[j+1:]
    return s1,s2

def CS_syms(name):
    return var('C{0},S{0}'.format(name))

def sym_less(A,B):
    return A.count_ops() + A.count(Symbol) < B.count_ops() + B.count(Symbol)

def get_angles(expr):
    angles = set()
    for s in expr.atoms(sin,cos):
        angles |= s.atoms(Symbol)
    return angles

from sympy import sin, cos, Symbol
from itertools import product
class Symoro:
    """Symbol manager, responsible for symbol replacing, file writing."""
    def __init__(self, file_out = 'disp'):
        """Default values correspond to empty dictionary and screen output."""

        self.file_out = file_out
        """Output descriptor. Can be None, 'disp', file
        defines the output destination"""

        self.sydi = {}
        """Dictionary. All the substitutions are saved in it"""

    def C4S4_simp(self,old_sym):
        """
        Example
        =======
        >> print C2S2_simp(sympify("-C**2*RL + S*(D - RL*S)"))
        D*S - RL
        """
        if old_sym.is_constant():
            return old_sym
        Res = old_sym.expand()
        for name in get_trig_pow_names(Res,4):
            C,S = CS_syms(name)
            E = Res
            EC = get_max_coef(E, C**4)
            ES = get_max_coef(E, S**4)
            if sym_less(EC,ES) or (EC/ES).is_constant() and EC/ES < 1:
                ES = EC
            E += ES - (ES*2*S**2*C**2).expand() -  (ES*S**4 + ES*C**4).expand()
            if sym_less(E,Res):
                Res = E
        return Res

    def C2S2_simp(self,old_sym):
        """
        Example
        =======
        >> print C2S2_simp(sympify("-C**2*RL + S*(D - RL*S)"))
        D*S - RL
        """
        if old_sym.is_constant():
            return old_sym
        Res = old_sym.expand()
        for name in get_trig_pow_names(Res):
            C,S = CS_syms(name)
            E = Res
            EC = get_max_coef(E, C**2)
            ES = get_max_coef(E, S**2)
            if sym_less(EC,ES) or (EC/ES).is_constant() and EC/ES < 1:
                ES = EC
            E += ES -  (ES * S**2 + ES * C**2).expand()
            if sym_less(E,Res):
                Res = E
        return Res

    def CS12_simp(self,old_sym, short_form = True):
        if old_sym.is_constant():
            return old_sym
        if short_form:
            c_names, s_names = get_trig_couple_names(old_sym)
            names = c_names & s_names
        else:
            names = get_angles(old_sym)
        orig = old_sym.expand()
        Res = orig
        for n1, n2 in product(names,names):
            if short_form:
                C1,S1 = CS_syms(n1)
                C2,S2 = CS_syms(n2)
                np1,nm1 = get_pos_neg(n1)
                np2,nm2 = get_pos_neg(n2)
                n12 = ang_sum(np1,np2,nm1,nm2)
                nm12 = ang_sum(np1,nm2,nm1,np2)
                C12,S12 = CS_syms(n12)
                C1m2,S1m2 = CS_syms(nm12)
            else:
                C1,S1 = cos(n1), sin(n1)
                C2,S2 = cos(n2), sin(n2)
                C12,S12 = cos(n1+n2), sin(n1+n2)
                C1m2,S1m2 = cos(n1-n2), sin(n1-n2)

            def try_opt(A,B,C,Res,E,sym_prev = None):
                EB = get_max_coef(E,B)
                EC = get_max_coef(E,C)
                if EB != 0 and EC != 0:
                    if sym_less(EC,EB):
                        EB = EC
                    E += (EB*A) -  (EB * (B + C)).expand()
                    if sym_less(E, Res):
                        Res = E
                        sym_prev = (A, B + C)
                return Res,sym_prev
            Res,sym = try_opt(S12, S1*C2, C1*S2, Res, orig)
            Res,sym = try_opt(S1m2, S1*C2, -C1*S2, Res, orig, sym)
            if sym_less(Res, orig):
                orig = Res
                if short_form:
                    self.add_to_dict(sym[0],sym[1])
            Res,sym = try_opt(C12, C1*C2, -S1*S2, Res, orig)
            Res,sym = try_opt(C1m2, C1*C2, S1*S2, Res, orig, sym)
            if sym_less(Res, orig):
                orig = Res
                if short_form:
                    self.add_to_dict(sym[0],sym[1])
        return Res

    def poly_simp(self,sym):
        """
        Tries to simplify A*C+A*D as A*(C+D)
        """
        sym_res = sym
        terms = sym.as_ordered_terms()
        for i in  range(1,len(terms)):
            for j in range(i):
                ratio = (terms[i]/terms[j]).as_numer_denom()
                if ratio[0].count(Symbol)  == terms[i].count(Symbol):
                    continue
                C = ratio[0]
                D = ratio[1]
                A = terms[i]/C
                if A.count(-1) == 1:
                    A = -A
                    D = -D
                    C = -C
                A = self.poly_simp(A)
                CD = C + D
                sym_new = sym - (A*C+A*D) + A*CD
                if sym_less(sym_new, sym_res):
                    sym_res = sym_new
        if sym_less(sym_res, sym):
            sym_res = self.poly_simp(sym_res)
        return sym_res

    def apply_mat(self,f,M,*args):
        """ Applies fucntion f to each element of M
        """
        for i in range(M.shape[0]):
            for j in range(M.shape[1]):
                if not M[i,j].is_constant():
                    if len(args) == 0:
                        M[i,j] = f(M[i,j])
                    else:
                        M[i,j] = f(M[i,j],*args)

    def add_to_dict(self,new_sym,old_sym):
        """Internal function.
        Extends symbol dictionary by (new_sym,old_sym) pair
        """
        if new_sym not in self.sydi.keys():
            self.sydi[new_sym] = old_sym
            self.write_equation(new_sym,old_sym)

    def trig_replace(self,M,angle,name):
        """Replaces trigonometric expressions cos(x)
        and sin(x) by CX and SX

        Parameters
        ==========
        M : var or Matrix
            Object of substitution
        angle : var
            symbol that stands for the angle value
        name : int or string
            brief name X for the angle

        Notes
        =====
        The cos(x) and sin(x) will be replaced by CX and SX,
        where X is the name and x is the angle
        """
        if sympify(angle).is_constant():
            return M
        cos_sym, sin_sym = CS_syms(name)
        sym_list = [(cos_sym,cos(angle)),(sin_sym,sin(angle))]
        subs_dict = {}
        for sym,sym_old in sym_list:
            subs_dict[sym_old] = sym
            if not sym in self.sydi:
                self.add_to_dict(sym, sym_old)
        for i1 in range(M.shape[0]):
            for i2 in range(M.shape[1]):
                M[i1,i2] = M[i1,i2].subs(subs_dict)
        return M

    def replace(self,old_sym,name,index = '',forced = False):
        """Creates a new symbol for the symbolic expression old_sym.

        Parameters
        ==========
        old_sym : var
            Symbolic expression to be substituted
        name : string
            denotion of the expression
        index : int or string, optional
            will be attached to the name. Usualy used for link or joint number.
            Parameter exists for usage convenience
        forced : bool, optional
            If True, the new symbol will be created even if old symbol
            is a simple expression

        Notes
        =====
        Generaly only complex expressions, which contain + - * / ** operations
        will be replaced by a new symbol
        """
        is_simple = old_sym.count_ops() == 0 or (-old_sym).count_ops() == 0
        if is_simple and not forced:
            return old_sym
        for i in (1, - 1):
            if not is_simple and i*old_sym in self.sydi.values():
                old_sym = i*self.sydi.keys()[self.sydi.values().index(i*old_sym)]
                if not forced:
                    return old_sym
                break
        new_sym = var(name + str(index))
        self.add_to_dict(new_sym, old_sym)
        return new_sym

    def mat_replace(self,M,name,index = '',
                    forced = False,
                    skip = 0,
                    symmet = False):
        """Replaces each element in M by symbol

        Parameters
        ==========
        M : Matrix
            Object of substitution
        name : string
            denotion of the expression
        index : int or string, optional
            will be attached to the name. Usualy used for link or joint number.
            Parameter exists for usage convenience
        forced : bool, optional
            If True, the new symbol will be created even if old symbol
            is a simple expression
        skip : int, optional
            Number of bottom rows of the matrix, which will be skipped.
            Used in case of Transformation matrix and forced = True.
        symmet : bool, optional
            If true, only for upper triangle part of the matrix
            symbols will be created. The bottom triangle part the
            same symbols will be used


        Returns
        =======
        M : Matrix
            Matrix with all the elements replaced

        Notes
        =====
        -Each element M_ij will be replaced by symbol name + i + j + index
        -There are two ways to use this function (examples):
            1)  >>> A = B+C+...
                >>> symo.mat_replace(A,'A')
                #for the case when expression B+C+... is too big
            2)  >>> A = symo.mat_replace(B+C+...,'A')
                #for the case when B+C+... is small enough
        """
        for i1 in range(M.shape[0] - skip):
            for i2 in range(M.shape[1]):
                if symmet and i2 < i1:
                    M[i1,i2] = M[i2,i1]
                    continue
                if M.shape[1] > 1:
                    name_index = name + str(i1 + 1) + str(i2 + 1)
                else:
                    name_index = name + str(i1 + 1)
                M[i1,i2] = self.replace(M[i1,i2],name_index,index,forced)
        return M

    def unfold(self,expr):
        """Unfold the expression using the dictionary.

        Parameters
        ==========
        expr : symbolic expression
            Symbolic expression to be unfolded

        Returns
        =======
        expr : symbolic expression
            Unfolded expression
        """
        while self.sydi.keys() & expr.atoms():
            expr = expr.subs(self.sydi)
        return expr

    def write_param(self,name,header,param,N):
        """Low-level function for writing the parameters table

        Parameters
        ==========
        name : string
            the name of the table
        header : list
            the table header
        param : callable (int) : list
            returns the list of parameters for
            the particular row of the table
        N : int
            number of lines in the table
        """
        self.write_line(name)
        self.write_line(l2str(header))
        for j in range(N):
            self.write_line(l2str(param(j)))
        self.write_line()

    def write_geom_param(self,robo,title = ''):
        """Writes the geometric parameters table

        Parameters
        ==========
        robo : Robot
            Instance of the parameter container
        title : string, optional
            The document title. Not used in case of internal using
        """
        if title != '':
            self.write_line(title)
            self.write_line()
        self.write_param('Geometric parameters',robo.get_geom_head(),
                         robo.get_geom_param,robo.NF)

    def write_inert_param(self,robo, name = 'Dynamic inertia parameters'):
        """Writes the inertia parameters table

        Parameters
        ==========
        robo : Robot
            Instance of the parameter container
        name : string, optional
            The table name. Not used in case of internal using.
        """
        self.write_param(name,robo.get_inert_head(),
                 robo.get_inert_param,robo.NL)

    def write_dynam_param(self,robo,title):
        """Writes the geometric parameters table

        Parameters
        ==========
        robo : Robot
            Instance of the parameter container.
        title : string
            The document title.

        Notes
        =====
        The synamic model generation program can be started with this function
        """
        self.write_line(title)
        self.write_line()
        self.write_geom_param(robo)
        self.write_inert_param(robo)
        self.write_param('External forces and joint parameters',
                         robo.get_ext_dynam_head(),
                         robo.get_ext_dynam_param,robo.NL)
        self.write_param('Base velicities parameters',robo.get_base_vel_head(),
                 robo.get_base_vel,3)


    def write_equation(self, A, B):
        """Writes the equation A = B into the output

        Parameters
        ==========
        A : expression or var
            left-hand side of the equation.
        B : expression or var
            right-hand side of the equation
        """
        self.write_line(str(A) + ' = ' + str(B))

    def write_line(self, line = ''):
        """Writes string data into tha output with new line symbol

        Parameters
        ==========
        line : string, optional
            Data to be written. If empty, it adds an empty line
        """
        if self.file_out == 'disp':
            print line
        elif self.file_out != None:
            self.file_out.write(str(line) + '\n')

    def file_open(self, robo, ext):
        """
        Initialize file stream

        Parameters
        ==========
        robo : Robot instance
            provides the robot's name
        ext : string
            provides the file name extention
        """
        self.file_out = open('models\\' + robo.name + '_' + ext + '.txt','w')

    def topological_sort(self, var_list):
        no_incom_edges = set(var_list)
        edge_dict = dict((v,0) for v in var_list)
#        print 'build graph',var_list
        while len(var_list) > 0:
            v = var_list.pop()
            for atom in self.sydi[v].atoms(Symbol):
                no_incom_edges -= {atom}
                if atom not in edge_dict:
                    edge_dict[atom] = 1
                    if atom in self.sydi:
                        var_list.append(atom)
                else:
                    edge_dict[atom] += 1
#        print 'enumerate graph',no_incom_edges
        var_fringe = list(no_incom_edges)
        sorted_vars = []
        while len(var_fringe) > 0:
            v = var_fringe.pop()
            sorted_vars.append(v)
            for atom in self.sydi[v].atoms(Symbol):
                edge_dict[atom] -= 1
                if edge_dict[atom] <= 0:
                    var_fringe.append(atom)
                    if atom not in self.sydi:
#                        print 'new symbol', atom
                        self.sydi[atom] = sympify(1.)
        return sorted_vars

    def as_function(self,name,var_type,var_size,*args):
        _intend = '    '
        fun = 'from numpy import sin,cos,sign,array\n'
        fun += 'def ' + name + '_func(*args):\n'
        for i,var_list in enumerate(args):
            fun += _intend
            for v in var_list:
                v_str = str(v)
                v_str = v_str.replace('[0]','[nul]')
                v_str = v_str.replace(']\n[',',')
                fun += v_str
                if v != var_list[-1]:
                    fun += ','
            fun += ' = args[{0}]\n'.format(i)
        to_compute_stack = []

        def write_vec(ret,name,length):
            for j in range(length):
                v = sympify(name  + str(j+1))
                if self.sydi.get(v) == None:
                    ret += '0'
                elif sympify(self.sydi[v]).is_constant():
                    ret += str(self.sydi[v])
                else:
                    ret += str(v)
                    to_compute_stack.append(v)
                if j < length-1:
                    ret += ','
            return ret
        #generates return statement
        ret = 'return array(['
        if var_type == 'matrix':
            for i in range(var_size[0]):
                ret += '['
                ret = write_vec(ret,name+ str(i+1),var_size[1])
                ret += ']'
                if i < var_size[0]-1:
                    ret += ','
        elif var_type == 'vector':
            ret = write_vec(ret,name,var_size)
        ret += '])'

        to_compute = self.topological_sort(to_compute_stack)
        for v in reversed(to_compute):
            if re.search(r'[\s,\[]'+str(v),fun) == None:
                fun += _intend + str(v) + ' = ' + str(self.sydi[v]) + '\n'
        fun += _intend + ret
        print fun
        exec fun in self.sydi
        return self.sydi[name + '_func']

#d = sympify("C2*C3*C4**2*C5**2*C6**4*D3**2*RL4*S5 + 2*C2*C3*C4**2*C5**2*C6**2*D3**2*RL4*S5*S6**2 + C2*C3*C4**2*C5**2*D3**2*RL4*S5*S6**4 + C2*C3*C4**2*C6**4*D3**2*RL4*S5**3 + 2*C2*C3*C4**2*C6**2*D3**2*RL4*S5**3*S6**2 + C2*C3*C4**2*D3**2*RL4*S5**3*S6**4 + C2*C3*C5**2*C6**4*D3**2*RL4*S4**2*S5 + 2*C2*C3*C5**2*C6**2*D3**2*RL4*S4**2*S5*S6**2 + C2*C3*C5**2*D3**2*RL4*S4**2*S5*S6**4 + C2*C3*C6**4*D3**2*RL4*S4**2*S5**3 + 2*C2*C3*C6**2*D3**2*RL4*S4**2*S5**3*S6**2 + C2*C3*D3**2*RL4*S4**2*S5**3*S6**4 - C3*C4**2*C5**2*C6**4*D3*RL4**2*S23*S5 - 2*C3*C4**2*C5**2*C6**2*D3*RL4**2*S23*S5*S6**2 - C3*C4**2*C5**2*D3*RL4**2*S23*S5*S6**4 - C3*C4**2*C6**4*D3*RL4**2*S23*S5**3 - 2*C3*C4**2*C6**2*D3*RL4**2*S23*S5**3*S6**2 - C3*C4**2*D3*RL4**2*S23*S5**3*S6**4 - C3*C5**2*C6**4*D3*RL4**2*S23*S4**2*S5 - 2*C3*C5**2*C6**2*D3*RL4**2*S23*S4**2*S5*S6**2 - C3*C5**2*D3*RL4**2*S23*S4**2*S5*S6**4 - C3*C6**4*D3*RL4**2*S23*S4**2*S5**3 - 2*C3*C6**2*D3*RL4**2*S23*S4**2*S5**3*S6**2 - C3*D3*RL4**2*S23*S4**2*S5**3*S6**4")
##d = sympify("C2**5*C3**3*C4**2*D3**2*RL4*S5 - C2**5*C3**3*C4**2*D3*RL4**2*S3*S5 + C2**5*C3**3*D3**2*RL4*S4**2*S5 - C2**5*C3**3*D3*RL4**2*S3*S4**2*S5 + C2**5*C3*C4**2*D3**2*RL4*S3**2*S5 - C2**5*C3*C4**2*D3*RL4**2*S3**3*S5 + C2**5*C3*D3**2*RL4*S3**2*S4**2*S5 - C2**5*C3*D3*RL4**2*S3**3*S4**2*S5 - C2**4*C3**4*C4**2*D3*RL4**2*S2*S5 - C2**4*C3**4*D3*RL4**2*S2*S4**2*S5 - C2**4*C3**2*C4**2*D3*RL4**2*S2*S3**2*S5 - C2**4*C3**2*D3*RL4**2*S2*S3**2*S4**2*S5 + 2*C2**3*C3**3*C4**2*D3**2*RL4*S2**2*S5 - 2*C2**3*C3**3*C4**2*D3*RL4**2*S2**2*S3*S5 + 2*C2**3*C3**3*D3**2*RL4*S2**2*S4**2*S5 - 2*C2**3*C3**3*D3*RL4**2*S2**2*S3*S4**2*S5 + 2*C2**3*C3*C4**2*D3**2*RL4*S2**2*S3**2*S5 - 2*C2**3*C3*C4**2*D3*RL4**2*S2**2*S3**3*S5 + 2*C2**3*C3*D3**2*RL4*S2**2*S3**2*S4**2*S5 - 2*C2**3*C3*D3*RL4**2*S2**2*S3**3*S4**2*S5 - 2*C2**2*C3**4*C4**2*D3*RL4**2*S2**3*S5 - 2*C2**2*C3**4*D3*RL4**2*S2**3*S4**2*S5 - 2*C2**2*C3**2*C4**2*D3*RL4**2*S2**3*S3**2*S5 - 2*C2**2*C3**2*D3*RL4**2*S2**3*S3**2*S4**2*S5 + C2*C3**3*C4**2*D3**2*RL4*S2**4*S5 - C2*C3**3*C4**2*D3*RL4**2*S2**4*S3*S5 + C2*C3**3*D3**2*RL4*S2**4*S4**2*S5 - C2*C3**3*D3*RL4**2*S2**4*S3*S4**2*S5 + C2*C3*C4**2*D3**2*RL4*S2**4*S3**2*S5 - C2*C3*C4**2*D3*RL4**2*S2**4*S3**3*S5 + C2*C3*D3**2*RL4*S2**4*S3**2*S4**2*S5 - C2*C3*D3*RL4**2*S2**4*S3**3*S4**2*S5 - C3**4*C4**2*D3*RL4**2*S2**5*S5 - C3**4*D3*RL4**2*S2**5*S4**2*S5 - C3**2*C4**2*D3*RL4**2*S2**5*S3**2*S5 - C3**2*D3*RL4**2*S2**5*S3**2*S4**2*S5")
#print 'det\n',d
#d = symo.C4S4_simp(d)
#print 'C4S4\n',d
#d = symo.C2S2_simp(d)
#print 'C2S2\n',d
#d = symo.CS12_simp(d)
#print 'CS12\n',d
#d = symo.poly_simp(d)
#print 'poly\n',d

#E = sympify("C3*D3*RL4*S5*(C5**2*(C2*D3*(C4**2*(C6**4 + S6**4) + C6**2*(2*C4**2*S6**2 + S4**2*(C6**2 + 2*S6**2))) - C4**2*RL4*S23*S6**4) + C6**4*(-C4**2 - S4**2)*(C5**2*RL4*S23 - S5**2*(C2*D3 - RL4*S23)) + S6**2*(-2*C6**2*(RL4*S23*(C4**2*C5**2 + C4**2*S5**2 + C5**2*S4**2) - S5**2*(C2*C4**2*D3 + S4**2*(C2*D3 - RL4*S23))) - S6**2*(C2*D3 - RL4*S23)*(-C4**2*S5**2 + S4**2*(-C5**2 - S5**2))))")
#print E
#print E.expand()
#E2 = Symoro().C4S4_simp(E)
#print E2
#E3 = Symoro().C2S2_simp(E2)
#print E3
#E4 = Symoro().poly_simp(E3)
#print E4
#print Symoro().C2S2_simp(sympify("C23**2 + S23**2"))
#print Symoro().CS12_simp(sympify("C2*D3*S3m78 - C2m7*D8*S3 - C3*D8*S2m7 - C3m78*D3*S2 + D2*S3"))
#print Symoro().CS12_simp(sympify("-C2m7*D8*S3 - C3*D8*S2m7 + D2*S3 - D3*S278m3"))