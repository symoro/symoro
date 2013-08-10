"""
This module of SYMORO package provides geometric models' computation.

The core symbolic library is sympy.

Needed modules : symoro.py

ECN - ARIA1 2013
"""

from sympy import Matrix, zeros, eye, sin, cos, sympify, trigsimp
from symoro import Symoro, Robot, hat

z_axis = Matrix([0,0,1])

def get_kRj(robo,antRj,k,chainj):
    T = eye(3)
    all_paral = True
    for i in chainj:
        if i > k:
            T = antRj[i]*T
        if antRj[i].col(2) != z_axis and robo.ant[i] != -1:
            all_paral = False
    return T, all_paral

def z_paral(T):
    return T.col(2) == z_axis

def find_r12(robo,chainj,antRj,j):
    r1 = robo.NL
    r2 = robo.NL
    rot12 = eye(3)
    orthog = False
    for i in reversed(chainj):
        if robo.sigma[i] == 0:
            rot12 *= antRj[i]
            if r1 == robo.NL:
                r1 = i
            elif r2 == robo.NL and rot12.col(2) != Matrix([0,0,1]):
                r2 = i
                if Matrix([0,0,1]).dot(rot12.col(2)) == 0:
                    orthog = True
                break
    return r1,r2,orthog

def get_sna(T):
    """Extracts the s,n,a vector basis of rotation 3x3 matrix
    from 4x4 transformation matrix

    Parameters
    ==========
    T : Matrix 4x4
        Transformation matrix

    Returns
    =======
    s : Matrix 3x1
    n : Matrix 3x1
    a : Matrix 3x1
    """
    R = get_r(T)
    return R.col(0), R.col(1), R.col(2)

def get_r(T):
    """Extracts rotation 3x3 matrix from 4x4 transformation matrix

    Parameters
    ==========
    T : Matrix 4x4
        Transformation matrix

    Returns
    =======
    get_r : Matrix 3x3
    """
    return T[:3,:3]

def get_p(T):
    """Extracts translation vector from 4x4 transformation matrix

    Parameters
    ==========
    T : Matrix 4x4
        Transformation matrix

    Returns
    =======
    get_p : Matrix 3x1
    """
    return T[:3,3]

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

def compute_transform(robo,symo,j,antRj,antPj):
    """Internal function. Computes rotation matrix and translation vector
    of ant_T_j homogenuous transform. Does the trigonometric subsctitution
    and saves the symbols into symo.sydi

    Notes
    =====
    antPj and antRj are the output parameters
    """
    antTj = transform(robo,j)
    for angle,name in robo.get_angles(j):
        antTj = symo.trig_replace(antTj,angle,name)
    antRj[j] = symo.mat_replace(get_r(antTj),'A',robo.num[j])
    antPj[j] = symo.mat_replace(get_p(antTj),'L',robo.num[j])

def compute_screw_transform(robo,symo,j,antRj,antPj,jTant):
    """Internal function. Computes the screw transformation matrix
    between ant[j] and j frames.

    Notes
    =====
    jTant is an output parameter
    """
    jRant = antRj[j].T
    ET = symo.mat_replace( - jRant*hat(antPj[j]),'JPR',robo.num[j])
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3,3).row_join(jRant)]))

def dgm_serial(robo, symo, k, j, fast_form = True, termin = False, all_trans = False):
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
        u.reverse()
        simp_further = lambda e2: robo.alpha[u[u.index(e2)+1]] == 0
    else:
        simp_further = lambda e2: robo.alpha[e2] == 0
    T_res = eye(4)
    T = eye(4)
    T_dict = {}
    for e2 in u:
        T = transform(robo,e2,invert)*T
        for ang,name in robo.get_angles(e2):
            symo.trig_replace(T,ang,name)
        T = T.applyfunc(symo.CS12_simp)
        if u[-1] != e2 and simp_further(e2):
            continue
        T_res = T*T_res
        T = eye(4)
        #make substitution with saving the var to the dictionary
        if fast_form:
            if robo.ant[e2] != k or not termin:
                name = 'U{0}T{1}'.format(robo.num[robo.ant[e2]],robo.num[j])
                symo.mat_replace(T_res,name)
        if all_trans:
            T_dict[(e2,j)] = T_res
    if all_trans: return T_dict
    else: return T_res

def dgm2(robo,symo,i,j,fast_form = True, forced = False):
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

    Notes
    =====
    Might be nesessary to make a symbol replacement after
    """
    if i == j:
        return eye(4)
    k = robo.common_root(i,j)
#    print i+1,j+1,k+1
    kTj = dgm_serial(robo,symo,k,j,fast_form,i==k)
    iTk = dgm_serial(robo,symo,i,k,fast_form,j==k)
    iTj = iTk*kTj
    paral = all(robo.alpha[indx] == 0
                for indx in robo.chain(i,k)[:-1] + robo.chain(j,k)[:-1])
    if k != i and k != j:
        ali_base = robo.alpha[robo.chain(i,k)[-1]]
        alj_base = robo.alpha[robo.chain(j,k)[-1]]
        base_paral = sin(ali_base - alj_base) == 0
    else:
        base_paral = True
    if paral and base_paral:
        iTj = iTj.applyfunc(symo.CS12_simp)
#        iTj = iTj.applyfunc(symo.poly_simp)
    if forced:
        name = 'T{0}T{1}'.format(robo.num[i],robo.num[j])
        symo.mat_replace(iTj,name,forced = forced)
    return iTj

def trans_name(robo,i,j,pattern = 'T{0}T{1}'):
    return 'T{0}T{1}'.format(robo.num[i],robo.num[j])

def dgm(robo, symo, i, j, fast_form = True, all_trans = False,
        forced = False, trig_subs = True):
    """Low-level Direct Geometric Model.

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
    k = robo.common_root(i,j)
    u = robo.chain(j,k)
    u2 = robo.chain(i,k)
    u.append(k)
    u.extend(reversed(u2))
    T_res = eye(4)
    T = eye(4)
    T_dict = {(j,j):T_res}
    for indx, e2 in enumerate(u):
        if u[-1] == e2:
            break
        e1 = u[indx+1]
        if e2 > e1:
            T = transform(robo,e2)*T
            if e1 != k or e1 == u[-1] :
                next_alpha = robo.alpha[e2]
            else:
                next_alpha = robo.alpha[e2] - robo.alpha[u[indx+2]]
        else:
            T = transform(robo,e1,True)*T
            if e1 != u[-1]:
                next_alpha = robo.alpha[u[indx+2]]
            else:
                next_alpha = 0
        if trig_subs:
            for ang,name in set(robo.get_angles(e2)) | set(robo.get_angles(e1)):
    #            print T,ang,name
                symo.trig_replace(T,ang,name)
            symo.apply_mat(symo.CS12_simp, T)
        else:
#                   symo.apply_mat(symo.CS12_simp, T,False)
            T = T.applyfunc(trigsimp)
        T_tmp = T*T_res
        T_dict[(e1,j)] = T_tmp
        if e1 != u[-1] and next_alpha == 0:
            continue
        T_res = T_tmp
        T = eye(4)
        #make substitution with saving the var to the dictionary
        if fast_form:
            if e1 != i:
                name = trans_name(robo,e1,j,'U{0}T{1}')
                symo.mat_replace(T_res,name)
            else:
                name = trans_name(robo,e1,j)
                symo.mat_replace(T_res,name,forced = forced)
    if all_trans: return T_dict
    else:
        return T_res

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

#TODO: validate for different structures
def direct_geometric(robo,i,j,fast_form = True):
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
    dgm(robo,symo, i - 1, j - 1, fast_form, write_res = True)
    symo.file_out.close()
    return symo.sydi

#robo = Robot.RX90()
#symo = Symoro()
##print 'Direct geometric model'
#T = dgm(robo,symo,-1,4,fast_form = False)
#print T
#print rot_trans('z', robo.theta[6], robo.r[6])* rot_trans('z', - robo.theta[1], - robo.r[1])
#print transform(robo,6,True),transform(robo,1),transform(robo,6,True)*transform(robo,1)
#T2 = dgm(robo,symo,8,9,fast_form = True)
#print T2

#robo = Robot.RX90()
#symo = Symoro(None)
#i,j = -1,3
#dgm(robo,symo,i,j, forced = True)
#
#f = symo.as_function(trans_name(robo,i,j),'matrix',[4,4],robo.q)
##dgm(robo,symo,j,i, forced = True)
##
##f2 = symo.as_function(trans_name(robo,j,i),'matrix',[4,4],robo.q)
#from numpy import dot
#print f([1,1,1,1,1,1])
#print f([1.5707963267948966, 1.5707963267948968, 1.0471975511965976, -1.5707963267948966, -1.0471975511965979, 3.141592653589793])
#print f([0,0,0,0,0,0])
