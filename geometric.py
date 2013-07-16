"""
This module of SYMORO package provides geometric models' computations.

The core symbolic library is sympy.

ECN - ARIA1 2013
"""

from sympy import Matrix,zeros,eye,sin,cos
from symoro import Symoro, RX90

def transform(robo,j,invert=False):
    """Transform matrix between frames j and ant[j]
    
    Parameters
    ==========
    j : int
        Frame index.
    invert : bool, optional
        Defines the transformation direction
        
    Returns
    =======
    transform : Matrix 4x4
        Transformation matrix. If invert is True then j_T_ant,
        else ant_T_j.
    """
    if not invert:
        R1 = rot_trans('z',robo.gamma[j],robo.b[j])
        R2 = rot_trans('x',robo.alpha[j],robo.d[j])
        R3 = rot_trans('z',robo.theta[j],robo.r[j])
        return R1*R2*R3
    else:
        R1 = rot_trans('z', - robo.gamma[j], - robo.b[j])
        R2 = rot_trans('x', - robo.alpha[j], - robo.d[j])
        R3 = rot_trans('z', - robo.theta[j], - robo.r[j])
        return R3*R2*R1    

def dgm_serial(robo,symo,k,j,fast_form = True):
    """Low-level Direct Geometric Model. For serial structures only.
    Used in Robot.dgm.
    
    Parameters
    ==========
    symo : Symoro
        Instance of Symoro. All the substitutions will be put into symo.sydi
    k : int
        To-frame index.
    j : int
        From-frame index.
    fast_form : bool, optional
        if False, result will be in unfolded mode (triginimetric 
        substitutions only)
        
    Returns
    =======
    T_res : Matrix 4x4
        Transformation matrix k_T_j
    """
    invert = (k > j)
    u = robo.chain(max(j,k),min(j,k))
    if invert:
        u = reversed(u)
    T_res = eye(4)
    T = eye(4)
    syms_CS = []
    names = []
    for e2 in u:
        Te = transform(robo,e2,invert)
        if invert:
            T = T*Te
        else:
            T = Te*T
        for ang,name in robo.get_angles(e2):
            T = symo.trig_replace(T,ang,name,syms = syms_CS)
            names.append(name)
            if len(names) == 2 and len(syms_CS) == 4:
                name = names[1] + names[0]
                syms_CS =symo.mat_CS12_simp(T,syms_CS,name)
                names = [name]
        if robo.alpha[e2] == 0 and robo.ant[e2] != k:
            continue           
        if invert:
            T_res = T_res*T
        else:
            T_res = T*T_res
        T = eye(4)
        syms_CS = []
        names = []
        #make substitution with saving the var to the dictionary
        if robo.ant[e2] != k and fast_form: 
            name = 'U{0}T{1}'.format(robo.num[robo.ant[e2]],robo.num[j])
            symo.mat_replace(T_res,name)
    name = 'T{0}T{1}'.format(robo.num[k],robo.num[j])
    symo.mat_replace(T_res,name,forced = True,skip = 1)
    return T_res

def dgm(robo,symo,i,j,fast_form = True):
    """Direct Geometric Model.
    
    Parameters
    ==========
    symo : Symoro
        Instance of Symoro. All the substitutions will be put into symo.sydi
    i : int
        To-frame index.
    j : int
        From-frame index.
    fast_form : bool, optional
        if False, result will be in unfolded mode (triginimetric 
        substitutions only)
        
    Returns
    =======
    T_res : Matrix 4x4
        Transformation matrix k_T_j
    """
    k = robo.common_root(i,j)
    kTj = dgm_serial(robo,symo,k,j,False,fast_form)
    iTk = dgm_serial(robo,symo,k,i,True,fast_form)
    return iTk*kTj
    
#section: GEOMETRIC
def rot(axis = 'z',th = 0):
    """Rotation matrix about axis
    
    Parameters
    ==========
    axis : {'x','y','z'}
        Rotation axis
    th : var
        Rotation angle
        
    Returns
    =======
    rot : Matrix 3x3
    """
    if axis == 'x':
        return  Matrix([[1,0,0],
                        [0,cos(th), - sin(th)],
                        [0,sin(th),cos(th)]])
    elif axis == 'y':
        return  Matrix([[cos(th),0,sin(th)],
                        [0,1,0],
                        [ - sin(th),0,cos(th)]])
    else:
        return  Matrix([[cos(th), - sin(th),0],
                        [sin(th),cos(th),0],
                        [0,0,1]])

def trans_vect(axis = 'z',p = 0):
    """Translation vector along axis
    
    Parameters
    ==========
    axis : {'x','y','z'}
        Translation axis
    p : var
        Translation distance
        
    Returns
    =======
    v : Matrix 3x1
    """
    axis_dict = {'x':0,'y':1,'z':2}
    v = zeros(3,1)
    v[axis_dict[axis]] = p
    return v

def trans(axis = 'z',p = 0):
    """Translation matrix along axis
    
    Parameters
    ==========
    axis : {'x','y','z'}
        Translation axis
    p : var
        Translation distance
        
    Returns
    =======
    trans : Matrix 4x4
    """
    return Matrix([eye(3).row_join(trans_vect(axis,p)),
                       [0,0,0,1]])

def rot_trans(axis = 'z',th = 0,p = 0):
    """Transformation matrix with rotation about and translation along axis
    
    Parameters
    ==========
    axis : {'x','y','z'}
        Transformation axis
    p : var
        Translation distance
    th : var
        Rotation angle
        
    Returns
    =======
    rot_trans : Matrix 4x4
    """
    return Matrix([rot(axis,th).row_join(trans_vect(axis,p)),
                       [0,0,0,1]])

def direct_geometric(robo,i,j,fast = True):
    """Computes trensformation matrix iTj.
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    i : int 
        the to-frame
    j : int
        the from-frame
    fast : bool
        if false, then the expressions will be unfolded
        
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    symo = Symoro()    
    symo.file_open(robo,'dgm')
    symo.write_geom_param(robo,'Direct Geometrix model')
    dgm_serial(robo,symo, i - 1, j - 1, fast)
    symo.file_out.close()
    return symo.sydi

robo = RX90()

print 'Direct geometric model'
direct_geometric(robo,0,6)
