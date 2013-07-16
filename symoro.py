"""
This module of SYMORO package provides description of the robot parametrizaion 
container and symbol replacer class.

The core symbolic library is sympy.

ECN - ARIA1 2013
"""
from sympy import Matrix,zeros,sin,cos,var,sympify,pi,sign,Integer

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
    sigma = None    
    """  joint type : list of int"""
    ant = None  
    """  index of antecedent joint : list of int"""
    num = None 
    """order numbers of joints (for display purposes) : list of int.
                The last number corresponds to the base frame -1 """
    
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
                    
def CartPole():
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
    
  
def RX90():
    """Generates Robot instance of RX90"""
    robo = Robot()
    #table of geometric parameters RX90
    robo.name = 'RX90'
    robo.NJ = 6
    robo.NL = 6
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
    robo.q = [var('Q{0}'.format(i)) for i in robo.num]
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
    robo.GAM = [0 for i in robo.num]
    robo.J = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                        'XY{0},YY{0},YZ{0},'
                        'XZ{0},YZ{0},ZZ{0}').format(i))) for i in robo.num]
    robo.G = Matrix([0,0,var('G3')])
    robo.num.append(0)  
    return robo

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
    
class Symoro:  
    """Symbol manager, responsible for symbol replacing, file writing."""
    def __init__(self):
        """Default values correspond to empty dictionary and screen output."""
        
        self.file_out = 'disp'
        """Output descriptor. Can be None, 'disp', file
        defines the output destination"""
        
        self.sydi = {} 
        """Dictionary. All the substitutions are saved in it"""
        
    def mat_C2S2_simp(self,M,name):
        """Siplify the trigonometric identity cos^2+sin^2=1
    
        Parameters
        ==========
        M : Matrix
            Object of simplification
        name : int or string
            brief angle name
            
        Notes
        =====
        In the expression to simplify sin(x) and cos(x) must be replaced
        with symbols SX and CX, where X is the name for x
        """        
        C,S = var('C{0},S{0}'.format(name))
        for i1 in range(M.shape[0]):
            for i2 in range(M.shape[1]):
                M[i1,i2] = M[i1,i2].subs(C**2 + S**2,Integer(1))
    
    def mat_CS12_simp(self,M,syms_CS,name):
        """Siplifies the trigonometric identities for cos(x+y) and sin(x+y)
    
        Parameters
        ==========
        M : Matrix
            Object of simplification
        syms_CS : list with pattern [C1, S1, C2, S2] 
            List of symbols of sin and cos of angles x and y
        name : int or string
            brief name for sum of angles
            
        Notes
        =====
        The cos(x+y) and sin(x+y) will be expressed through C1, S1, C2, S2
        what is computationaly more efficient.
        """
        C1,S1,C2,S2 = syms_CS
        C12,S12 = var('C{0},S{0}'.format(name))
        self.sydi[S12] = C1*S2 + S1*C2
        self.sydi[C12] = C1*C2 - S1*S2
        self.write_equation(C12,self.sydi[C12])
        self.write_equation(S12,self.sydi[S12])
        for i1 in range(M.shape[0]):
            for i2 in range(M.shape[1]):
                M[i1,i2] = M[i1,i2].subs(self.sydi[S12],S12)
                M[i1,i2] = M[i1,i2].subs(self.sydi[C12],C12)
        return [C12,S12]
    
    def trig_replace(self,M,angle,name,syms = None):
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
        syms : list, optional
            both CX and SX will be stored there. 
            Used for further trigonometric simplification.
            Order is important.
            
        Notes
        =====
        The cos(x) and sin(x) will be replaced by CX and SX,
        where X is the name and x is the angle
        """
        if sympify(angle).is_constant():
            return M
        cos_sym = var('C' + str(name))
        sin_sym = var('S' + str(name))
        sym_list = [(cos_sym,cos(angle)),(sin_sym,sin(angle))]
        subs_dict = {}
        for sym,sym_old in sym_list:
            subs_dict[sym_old] = sym
            if not sym in self.sydi:
                self.sydi[sym] = sym_old        
                self.write_equation(sym,self.sydi[sym])
            if syms != None:
                syms.append(sym)
        return M.subs(subs_dict)
    
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
        is_complex = old_sym.count_ops() == 0 or ( - old_sym).count_ops() == 0
        if is_complex and not forced:
            return old_sym        
        for i in (1, - 1):
            if i*old_sym in self.sydi.values():
                old_sym = i*self.sydi.keys()[self.sydi.values().index(i*old_sym)] 
                if not forced:
                    return old_sym
                break
        new_sym = var(name + str(index))
        self.sydi[new_sym] = old_sym
        self.write_equation(new_sym,old_sym)
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
            1)  A = B+C+...
                symo.mat_replace(A,'A')
                #for the case when expression B+C+... is too big
            2)  A = symo.mat_replace(B+C+...,'A')
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
        while any([self.sydi.has_key(a) for a in expr.atoms()]):
            expr = expr.subs(self.sydi)
        return expr
    
    def mat_unfold(self,M):
        """Unfold each elemet in the matrix.
    
        Parameters
        ==========
        M : Matrix
            Matrix of symbolic expressions to be unfolded
        
        See Also
        ========
        unfold
        """
        for i1 in range(M.shape[0]):
            for i2 in range(M.shape[1]):
                M[i1,i2] = self.unfold(M[i1,i2])
    
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
                         robo.get_geom_param,robo.NJ)
                         
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
                 robo.get_inert_param,robo.NJ)
    
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
                         robo.get_ext_dynam_param,robo.NJ)
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
    