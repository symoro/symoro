from sympy import Matrix,eye,zeros,sin,cos,var,trigsimp,sympify,pi,sign

axis_dict = {'x':0,'y':1,'z':2}

class Robot:
    NL = None
    NJ = None
    sigma = None
    ant = None
    num = None
    theta = None
    r = None
    alpha = None
    d = None
    gamma = None
    b = None
    J = None
    MS = None
    M = None
    G = None 
    GAM = None
    W0 = None
    WP0 = None
    V0 = None
    VP0 = None
    qdot = None
    qddot = None
    n_ext = None
    f_ext = None
    FS = None
    FV = None
    IA = None

    def fric_v(self,j):
        return self.FV[j]*self.qdot[j]
    
    def fric_s(self,j):
        return self.FS[j]*sign(self.qdot[j])

    def tau_ia(self,j):
        return self.IA[j]*self.qddot[j]
    
    def transform(self,index,invert=False):
        if not invert:
            R1 = rot_trans('z',self.gamma[index],self.b[index])
            R2 = rot_trans('x',self.alpha[index],self.d[index])
            R3 = rot_trans('z',self.theta[index],self.r[index])
            return R1*R2*R3
        else:
            R1 = rot_trans('z',-self.gamma[index],-self.b[index])
            R2 = rot_trans('x',-self.alpha[index],-self.d[index])
            R3 = rot_trans('z',-self.theta[index],-self.r[index])
            return R3*R2*R1
            
    def get_angles(self,index):
        return (self.gamma[index],self.alpha[index],self.theta[index])

def CartPole():
    robo = Robot()
    robo.ant = (-1,0)
    robo.sigma = (1,0)
    robo.alpha = (pi/2,pi/2)
    robo.d = (0,0)
    robo.theta = (pi/2,var('th'))
    robo.r = (var('x'),0)
    robo.b = (0,0)
    robo.gamma = (0,0)
    robo.num = range(1,3)
    robo.NJ = 2
    robo.NL = 2
    robo.n_ext = [zeros(3,1) for i in robo.num]
    robo.f_ext = [zeros(3,1) for i in robo.num]
    robo.FS = [0 for i in robo.num]
    robo.IA = [0 for i in robo.num]
    robo.FV = [0 for i in robo.num]
    robo.MS = [zeros(3,1) for i in robo.num]
    robo.MS[1][0] = var('m2l')
    robo.M = [var('m{0}'.format(i)) for i in robo.num]
    robo.GAM = [var('u'),0]
    robo.J = [zeros(3) for i in robo.num]
    robo.J[1][2,2] = var('m2ll')
    robo.G = Matrix([0,0,-var('g')])
    robo.W0 = zeros(3,1)
    robo.WP0 = zeros(3,1)
    robo.V0 = zeros(3,1)
    robo.VP0 = zeros(3,1)
    robo.q = var('x,th')
    robo.qdot =  var('xd,thd')
    robo.qddot = var('xdd,thdd')        
    return robo
    
  
def RX90():
    robo = Robot()
    #table of geometric parameters RX90
    robo.NJ = 6
    robo.NL = 6
    robo.num = range(1,robo.NJ+1)
    robo.ant = (-1,0,1,2,3,4)
    robo.sigma = (0,0,0,0,0,0)
    robo.alpha = (0,pi/2,0,-pi/2,pi/2,-pi/2)
    robo.d = (0,0,var('D3'),0,0,0)
    robo.theta = list(var('th1:7'))
    robo.r = (0,0,0,var('RL4'),0,0)
    robo.b = (0,0,0,0,0,0)
    robo.gamma = (0,0,0,0,0,0)    
    robo.W0 = zeros(3,1)
    robo.WP0 = zeros(3,1)
    robo.V0 = zeros(3,1)
    robo.VP0 = zeros(3,1)    
    robo.q = [var('Q{0}'.format(i)) for i in robo.num]
    robo.qdot =  [var('QP{0}'.format(i)) for i in robo.num]
    robo.qddot = [var('QDP{0}'.format(i)) for i in robo.num]    
    robo.n_ext= [zeros(3,1) for i in robo.num]
    robo.n_ext[-1] = Matrix(var('CX{0},CY{0},CZ{0}'.format(robo.num[-1])))
    robo.f_ext = [zeros(3,1) for i in robo.num]
    robo.f_ext[-1] = Matrix(var('FX{0},FY{0},FZ{0}'.format(robo.num[-1])))
    robo.FS = var('FS0:{0}'.format(robo.NL))
    robo.IA = var('IA0:{0}'.format(robo.NL))
    robo.FV = var('FV0:{0}'.format(robo.NL))
    robo.MS = [Matrix(var('MX{0},MY{0},MZ{0}'.format(i))) for i in robo.num]
    robo.M = [var('M{0}'.format(i)) for i in robo.num]
    robo.GAM = [0 for i in robo.num]
    robo.J = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                        'XY{0},YY{0},YZ{0},'
                        'XZ{0},YZ{0},ZZ{0}').format(i))) for i in robo.num]
    robo.G = Matrix([0,0,var('G3')])
    return robo

#section GEOMETRIC
def rot(axis = 'z',th = 0):
    if axis == 'x':
        return  Matrix([[1,0,0],
                        [0,cos(th),-sin(th)],
                        [0,sin(th),cos(th)]])
    elif axis == 'y':
        return  Matrix([[cos(th),0,sin(th)],
                        [0,1,0],
                        [-sin(th),0,cos(th)]])
    else:
        return  Matrix([[cos(th),-sin(th),0],
                        [sin(th),cos(th),0],
                        [0,0,1]])

def trans_vect(axis = 'z',p = 0):
    v = zeros(3,1)
    v[axis_dict[axis]] = p
    return v

def trans(axis = 'z',p = 0):
    return Matrix([eye(3).row_join(trans_vect(axis,p)),
                       [0,0,0,1]])

def rot_trans(axis = 'z',th = 0,p = 0):
    return Matrix([rot(axis,th).row_join(trans_vect(axis,p)),
                       [0,0,0,1]])
#endsection GEOMETRIC
                       
#section MISCELANOUS
def get_r(T):
    return T[:3,:3]

def get_p(T):
    return T[:3,3]

def hat(v):
    return Matrix([[0,-v[2],v[1]],
                   [v[2],0,-v[0]],
                   [-v[1],v[0],0]])

def mat_trigsimp(M):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            M[i1,i2] = trigsimp(M[i1,i2])

def trig_replace(M,sym_dict,angle,name,disp = True):
    if sympify(angle).is_constant():
        return M
    cos_sym = var('C'+str(name))
    sin_sym = var('S'+str(name))
    sym_dict[cos_sym] = cos(angle)
    sym_dict[sin_sym] = sin(angle)
    if disp:
        print cos_sym,'=',sym_dict[cos_sym]
        print sin_sym,'=',sym_dict[sin_sym]
    return M.subs({sym_dict[sin_sym]:sin_sym,sym_dict[cos_sym]:cos_sym})

def sym_replace(old_sym,sym_dict,name,index='',disp=True,forced=False):
    if old_sym.count_ops() != 0 and (-old_sym).count_ops() != 0 or forced:
        name_index = name+str(index)
        new_sym = var(name_index)
        sym_dict[new_sym] = old_sym
        if disp:
            print new_sym, '=', old_sym
        return new_sym
    else:
        return old_sym

def sym_mat_replace(M,sym_dict,name,index = '',disp = True,forced = False):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            if M.shape[1] > 1:
                name_index = name+str(i1+1)+str(i2+1)
            else:
                name_index = name+str(i1+1)
            M[i1,i2] = sym_replace(M[i1,i2],sym_dict,name_index,index,disp,forced)
    return M
    
def sym_mats_replace(M,sym_dict,name,index = '',disp = True,forced = False):
    for i2 in range(M.shape[1]):
        for i1 in range(i2,M.shape[0]):
            name_index = name+str(i1+1)+str(i2+1)
            M[i1,i2] = sym_replace(M[i1,i2],sym_dict,name_index,index,disp,forced)
            M[i2,i1] = M[i1,i2]
    return M

def unfold(expr,symb_dict):
    while any([symb_dict.has_key(a) for a in expr.atoms()]):
        expr = expr.subs(symb_dict)
    return expr

def mat_unfold(M,symb_dict):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            M[i1,i2] = unfold(M[i1,i2],symb_dict)
