"""
This module of SYMORO package provides geometric models' computation.

The core symbolic library is sympy.

Needed modules : symoro.py

ECN - ARIA1 2013
"""

from sympy import Matrix, zeros, eye, sin, cos
from symoro import Symoro, Robot, Init, hat

Z_AXIS = Matrix([0, 0, 1])

class Transform():
    @classmethod
    def sna(self, T):
        """Extracts the s, n, a vector basis of rotation 3x3 matrix
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
        R = Transform.R(T)
        return R.col(0), R.col(1), R.col(2)

    @classmethod
    def R(self, T):
        """Extracts rotation 3x3 matrix from 4x4 transformation matrix

        Parameters
        ==========
        T : Matrix 4x4
            Transformation matrix

        Returns
        =======
        get_r : Matrix 3x3
        """
        return T[:3, :3]

    @classmethod
    def P(self, T):
        """Extracts translation vector from 4x4 transformation matrix

        Parameters
        ==========
        T : Matrix 4x4
            Transformation matrix

        Returns
        =======
        get_p : Matrix 3x1
        """
        return T[:3, 3]

    @classmethod
    def kRj(self, robo, antRj, k, chainj):
        T = eye(3)
        all_paral = True
        for i in chainj:
            if i > k:
                T = antRj[i]*T
            if antRj[i].col(2) != Z_AXIS and robo.ant[i] != -1:
                all_paral = False
        return T, all_paral

    @classmethod
    def find_r12(self, robo, chainj, antRj, j):
        r1 = robo.NL
        r2 = robo.NL
        rot12 = eye(3)
        orthog = False
        for i in reversed(chainj):
            if robo.sigma[i] == 0:
                rot12 *= antRj[i]
                if r1 == robo.NL:
                    r1 = i
                elif r2 == robo.NL and rot12.col(2) != Matrix([0, 0, 1]):
                    r2 = i
                    if Matrix([0, 0, 1]).dot(rot12.col(2)) == 0:
                        orthog = True
                    break
        return r1, r2, orthog

    @classmethod
    def z_paral(self, T):
        return T.col(2) == Z_AXIS


def transform(robo, j, invert=False):
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
        R1 = rot_trans('z', robo.gamma[j], robo.b[j])
        R2 = rot_trans('x', robo.alpha[j], robo.d[j])
        R3 = rot_trans('z', robo.theta[j], robo.r[j])
        return R1*R2*R3
    else:
        R1 = rot_trans('z', -robo.gamma[j], -robo.b[j])
        R2 = rot_trans('x', -robo.alpha[j], -robo.d[j])
        R3 = rot_trans('z', -robo.theta[j], -robo.r[j])
        return R3*R2*R1

def compute_transform(robo, symo, j, antRj, antPj):
    """Internal function. Computes rotation matrix and translation vector
    of ant_T_j homogenuous transform. Does the trigonometric subsctitution
    and saves the symbols into symo.sydi

    Notes
    =====
    antPj and antRj are the output parameters
    """
    antTj = transform(robo, j)
    for angle, name in robo.get_angles(j):
        antTj = symo.trig_replace(antTj, angle, name)
    antRj[j] = symo.mat_replace(Transform.R(antTj), 'A', robo.num[j])
    antPj[j] = symo.mat_replace(Transform.P(antTj), 'L', robo.num[j])

def compute_screw_transform(robo, symo, j, antRj, antPj, jTant):
    """Internal function. Computes the screw transformation matrix
    between ant[j] and j frames.

    Notes
    =====
    jTant is an output parameter
    """
    jRant = antRj[j].T
    ET = symo.mat_replace( -jRant*hat(antPj[j]), 'JPR', robo.num[j])
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3, 3).row_join(jRant)]))

def trans_name(robo, i, j, pattern='T{0}T{1}'):
    return 'T{0}T{1}'.format(robo.num[i], robo.num[j])

def dgm_left(robo, symo, i, j, trig_subs=True):
    k = robo.common_root(i, j)
    chain1 = robo.chain(j, k)
    chain2 = robo.chain(i, k)
    chain2.reverse()
    complete_chain = (chain1 + chain2 + [None])
    T_out = {(j,j):eye(4)}
    T_res = eye(4)
    T = eye(4)
    for indx, x in enumerate(complete_chain[:-1]):
        inverted = indx >= len(chain1)
        T = transform(robo, x, inverted) * T
        if trig_subs:
            for ang, name in robo.get_angles(x):
                symo.trig_replace(T, ang, name)
        T = T.expand()
        T = T.applyfunc(symo.CS12_simp)
        x_next = complete_chain[indx + 1]
        if inverted: trans_name = (x,j)
        else: trans_name = (robo.ant[x],j)
        T_out[trans_name] = T * T_res
        if robo.paral(x, x_next):
            continue
        T_res = T * T_res
        T = eye(4)
    return T_out

def dgm_right(robo, symo, i, j, fast_form=True,
            forced=False, trig_subs=True):
    k = robo.common_root(i, j)
    chain1 = robo.chain(i, k)
    chain2 = robo.chain(j, k)
    chain2.reverse()
    complete_chain = (chain1 + chain2 + [None])
    T_out = {(i,i):eye(4)}
    T_res = eye(4)
    T = eye(4)
    for indx, x in enumerate(complete_chain[:-1]):
        inverted = indx >= len(chain1)
        T = T * transform(robo, x, inverted)
        if trig_subs:
            for ang, name in robo.get_angles(x):
                symo.trig_replace(T, ang, name)
        T = T.expand()
        T = T.applyfunc(symo.CS12_simp)
        x_next = complete_chain[indx + 1]
        if inverted: trans_name = (i,x)
        else: trans_name = (i,robo.ant[x])
        T_out[trans_name] = T_res * T
        if robo.paral(x, x_next):
            continue
        T_res = T_res * T
        T = eye(4)
    return T_out

def dgm_one(robo, symo, i, j, fast_form=True,
            forced=False, trig_subs=True):
    k = robo.common_root(i, j)
    chain1 = robo.chain(j, k)
    chain2 = robo.chain(i, k)
    chain2.reverse()
    complete_chain = (chain1 + chain2 + [None])
    T_res = eye(4)
    T = eye(4)
    for indx, x in enumerate(complete_chain[:-1]):
        inverted = indx >= len(chain1)
        T = transform(robo, x, inverted) * T
        if trig_subs:
            for ang, name in robo.get_angles(x):
                symo.trig_replace(T, ang, name)
        T = T.expand()
        T = T.applyfunc(symo.CS12_simp)
        x_next = complete_chain[indx + 1]
        if robo.paral(x, x_next):
            continue
        T_res = T * T_res
        T = eye(4)
        if fast_form:
            dgm_rename(robo, symo, T_res, x, i, j, inverted, forced)
    return T_res

def dgm_rename(robo, symo, T_res, x, i, j, inverted, forced):
    if inverted:
        name = trans_name(robo, x, j)
        forced_now = x == i
    else:
        name = trans_name(robo, robo.ant[x], j)
        forced_now = robo.ant[x] == i
    symo.mat_replace(T_res, name, forced=forced and forced_now)

def dgm(robo, symo, i, j, key='one', fast_form=True,
              forced=False, trig_subs=True):
    """must be the final DGM function

     Parameters
    ==========
    symo : Symoro
        Instance of Symoro. All the substitutions will
        be put into symo.sydi
    i : int
        To-frame index.
    j : int
        From-frame index.
    fast_form : bool, optional
        if False, result will be in unfolded mode (triginimetric
        substitutions only)
    """
    if key == 'left':
        return dgm_left(robo, symo, i, j, trig_subs)
    elif key == 'right':
        return dgm_right(robo, symo, i, j, trig_subs)
    else:
        return dgm_one(robo, symo, i, j, fast_form, forced, trig_subs)


#def dgm_old(robo, symo, i, j, fast_form=True, all_trans=False,
#        forced=False, trig_subs=True):
#    """Low-level Direct Geometric Model.
#
#    Parameters
#    ==========
#    symo : Symoro
#        Instance of Symoro. All the substitutions will
#        be put into symo.sydi
#    i : int
#        To-frame index.
#    j : int
#        From-frame index.
#    fast_form : bool, optional
#        if False, result will be in unfolded mode (triginimetric
#        substitutions only)
#
#    Returns
#    =======
#    T_res : Matrix 4x4
#        Transformation matrix k_T_j
#    """
#    k = robo.common_root(i, j)
#    chain1 = robo.chain(j, k)
#    chain2 = robo.chain(i, k)
#    chain2.reverse()
#    T_res = eye(4)
#    T = eye(4)
#    T_dict = {(j, j):T_res}
#    next_alpha = 1
#    for step_two, chain in enumerate((chain1,chain2)):
#        for indx, e2 in enumerate(chain):
#            T = transform(robo, e2, step_two)*T
#            if step_two and indx + 1 < len(chain):
#                next_alpha = robo.alpha[chain[indx+1]]
#            elif not step_two and e2 == chain[-1]:
#                if len(chain2):
#                    next_alpha = robo.alpha[e2] - robo.alpha[chain[-1]]
#                else:
#                    next_alpha = 1
#            elif step_two and e2 == chain[-1]:
#                next_alpha = 1
#            else:
#                next_alpha = robo.alpha[e2]
#            if trig_subs:
#                for ang, name in robo.get_angles(e2):
#                    symo.trig_replace(T, ang, name)
#            T = T.expand()
#            T = T.applyfunc(symo.CS12_simp)
#            T_tmp = T*T_res
#            if all_trans:
#                if step_two: T_dict[(e2, j)] = T_tmp
#                else: T_dict[(robo.ant[e2], j)] = T_tmp
#            if next_alpha == 0:
#                continue
#            T_res = T_tmp
#            T = eye(4)
#            if fast_form:
#                if step_two:
#                    name = trans_name(robo, e2, j)
#                    forced_now = e2 == i
#                else:
#                    name = trans_name(robo, robo.ant[e2], j)
#                    forced_now = robo.ant[e2] == i
#                symo.mat_replace(T_res, name, forced = forced and forced_now)
#    if all_trans: return T_dict
#    else: return T_res

def rot(axis='z', th=0):
    """Rotation matrix about axis

    Parameters
    ==========
    axis : {'x', 'y', 'z'}
        Rotation axis
    th : var
        Rotation angle

    Returns
    =======
    rot : Matrix 3x3
    """
    if axis == 'x':
        return  Matrix([[1, 0, 0],
                        [0, cos(th), - sin(th)],
                        [0, sin(th), cos(th)]])
    elif axis == 'y':
        return  Matrix([[cos(th), 0, sin(th)],
                        [0, 1, 0],
                        [ - sin(th), 0, cos(th)]])
    else:
        return  Matrix([[cos(th), - sin(th), 0],
                        [sin(th), cos(th), 0],
                        [0, 0, 1]])

def trans_vect(axis='z', p=0):
    """Translation vector along axis

    Parameters
    ==========
    axis : {'x', 'y', 'z'}
        Translation axis
    p : var
        Translation distance

    Returns
    =======
    v : Matrix 3x1
    """
    axis_dict = {'x':0, 'y':1, 'z':2}
    v = zeros(3, 1)
    v[axis_dict[axis]] = p
    return v

def trans(axis ='z', p=0):
    """Translation matrix along axis

    Parameters
    ==========
    axis : {'x', 'y', 'z'}
        Translation axis
    p : var
        Translation distance

    Returns
    =======
    trans : Matrix 4x4
    """
    return Matrix([eye(3).row_join(trans_vect(axis, p)),
                       [0, 0, 0, 1]])

def rot_trans(axis='z', th=0, p=0):
    """Transformation matrix with rotation about and
    translation along axis

    Parameters
    ==========
    axis : {'x', 'y', 'z'}
        Transformation axis
    p : var
        Translation distance
    th : var
        Rotation angle

    Returns
    =======
    rot_trans : Matrix 4x4
    """
    return Matrix([rot(axis, th).row_join(trans_vect(axis, p)),
                       [0, 0, 0, 1]])

def compute_rot_trans(robo, symo):
    #init transformation
    antRj = Init.init_mat(robo)
    antPj = Init.init_vec(robo)
    for j in xrange(robo.NL):
        compute_transform(robo, symo, j, antRj, antPj)
    return antRj, antPj

#TODO: validate for different structures
def direct_geometric(robo, i, j, fast_form=True):
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
    symo.file_open(robo, 'dgm')
    symo.write_geom_param(robo, 'Direct Geometrix model')
    dgm(robo, symo, i - 1, j - 1, fast_form, write_res=True)
    symo.file_out.close()
    return symo


#print 'Direct geometric model'
def a():
    robo = Robot.SR400()
    symo = Symoro()
    T = dgm(robo, symo, 5, 0, fast_form=True, trig_subs=True)
    return symo.gen_func('DGM_generated1', T)
def b():
    robo = Robot.SR400()
    symo = Symoro()
    T = dgm(robo, symo, 5, 0, fast_form=True, trig_subs=True)
    return symo.gen_func('DGM_generated2', T)

##from timeit import timeit
##print timeit(a, number = 10)
##print timeit(b, number = 10)
##
import profile
#
#profile.run('b()', sort = 'cumtime')
#profile.run('b()')
from timeit import timeit
from numpy import matrix

#print (a() + b()).expand()
f1 = a()
f2 = b()
print timeit(f1,number = 10000)
print timeit(f2,number = 10000)
#print matrix(f2([1,2,3,4,5,6]))*matrix(f1([1,2,3,4,5,6]))
print f1()*matrix(f2())
###print rot_trans('z', robo.theta[6], robo.r[6])* rot_trans('z', - robo.theta[1], - robo.r[1])
#print transform(robo, 6, True), transform(robo, 1), transform(robo, 6, True)*transform(robo, 1)
#T2 = dgm(robo, symo, 8, 9, fast_form = True)
#print T2
 #
#robo = Robot.RX90()
#symo = Symoro()
#i, j = -1, 5
#T = dgm_final(robo, symo, i, j)

#print robo.get_q_vec()

#f = symo.gen_func('DGM_generated', T, robo.get_q_vec())
#dgm(robo, symo, j, i, forced = True)
##
#f2 = symo.as_function(trans_name(robo, j, i), 'matrix', [4, 4], robo.get_q_vec())
#from numpy import dot
#print timeit(f,number = 1000)# 0.16
#print f([1.5707963267948966, 1.5707963267948968, 1.0471975511965976, -1.5707963267948966, -1.0471975511965979, 3.141592653589793])
#print f2([0, 0, 2, 0, 0, 0])*f([0, 0, 2, 0, 0, 0])

