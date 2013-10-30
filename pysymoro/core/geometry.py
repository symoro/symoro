"""
This module of SYMORO package provides geometric models' computation.

The core symbolic library is sympy.

Needed modules: symoro.py

ECN - ARIA1 2013
"""

from sympy import Matrix, zeros, eye, sin, cos
from pysymoro.core.symoro import Symoro, Init, hat

Z_AXIS = Matrix([0, 0, 1])


class Transform():
    @classmethod
    def sna(self, T):
        """Extracts the s, n, a vector basis of rotation 3x3 matrix
        from 4x4 transformation matrix

        Parameters
        ==========
        T: Matrix 4x4
            Transformation matrix

        Returns
        =======
        s: Matrix 3x1
        n: Matrix 3x1
        a: Matrix 3x1
        """
        R = Transform.R(T)
        return R.col(0), R.col(1), R.col(2)

    @classmethod
    def R(self, T):
        """Extracts rotation 3x3 matrix from 4x4 transformation matrix

        Parameters
        ==========
        T: Matrix 4x4
            Transformation matrix

        Returns
        =======
        get_r: Matrix 3x3
        """
        return T[:3, :3]

    @classmethod
    def P(self, T):
        """Extracts translation vector from 4x4 transformation matrix

        Parameters
        ==========
        T: Matrix 4x4
            Transformation matrix

        Returns
        =======
        get_p: Matrix 3x1
        """
        return T[:3, 3]

    @classmethod
    def kRj(self, robo, antRj, k, chainj):
        T = eye(3)
        all_paral = True
        for i in chainj:
            if i > k:
                T = antRj[i]*T
            if antRj[i].col(2) != Z_AXIS and robo.ant[i] != 0:
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


def _transform(robo, j, invert=False):
    """Transform matrix between frames j and ant[j]

    Parameters
    ==========
    j: int
        Frame index.
    invert: bool, optional
        Defines the transformation direction

    Returns
    =======
    transform: Matrix 4x4
        Transformation matrix. If invert is True then j_T_ant,
        else ant_T_j.
    """
    if not invert:
        R1 = _rot_trans('z', robo.gamma[j], robo.b[j])
        R2 = _rot_trans('x', robo.alpha[j], robo.d[j])
        R3 = _rot_trans('z', robo.theta[j], robo.r[j])
        return R1*R2*R3
    else:
        R1 = _rot_trans('z', -robo.gamma[j], -robo.b[j])
        R2 = _rot_trans('x', -robo.alpha[j], -robo.d[j])
        R3 = _rot_trans('z', -robo.theta[j], -robo.r[j])
        return R3*R2*R1


#TODO: rewrite the description
def _transform_const_sep(robo, j, invert=False):
    """Transform matrix between frames j and ant[j]

    Parameters
    ==========
    j: int
        Frame index.
    invert: bool, optional
        Defines the transformation direction

    Returns
    =======
    transform: Matrix 4x4
        Transformation matrix. If invert is True then j_T_ant,
        else ant_T_j.
    """
    if not invert:
        R1 = _rot_trans('z', robo.gamma[j], robo.b[j])
        R2 = _rot_trans('x', robo.alpha[j], robo.d[j])
        R3 = _rot_trans('z', th=robo.theta[j])
        R4 = _rot_trans('z', p=robo.r[j])
        return R1, R2, R3, R4
    else:
        R1 = _rot_trans('z', -robo.gamma[j], -robo.b[j])
        R2 = _rot_trans('x', -robo.alpha[j], -robo.d[j])
        R3 = _rot_trans('z', th=-robo.theta[j])
        R4 = _rot_trans('z', p=-robo.r[j])
        return R1, R2, R3, R4


def compute_transform(robo, symo, j, antRj, antPj):
    """Internal function. Computes rotation matrix and translation vector
    of ant_T_j homogenuous transform. Does the trigonometric subsctitution
    and saves the symbols into symo.sydi

    Notes
    =====
    antPj and antRj are the output parameters
    """
    antTj = _transform(robo, j)
    for angle, name in robo.get_angles(j):
        antTj = symo.trig_replace(antTj, angle, name)
    antRj[j] = symo.mat_replace(Transform.R(antTj), 'A', j)
    antPj[j] = symo.mat_replace(Transform.P(antTj), 'L', j)


def compute_screw_transform(robo, symo, j, antRj, antPj, jTant):
    """Internal function. Computes the screw transformation matrix
    between ant[j] and j frames.

    Notes
    =====
    jTant is an output parameter
    """
    jRant = antRj[j].T
    ET = symo.mat_replace(-jRant*hat(antPj[j]), 'JPR', j)
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3, 3).row_join(jRant)]))


def trans_name(robo, i, j, pattern='T{0}T{1}'):
    return 'T%sT%s' % (i, j)


def _dgm_left(robo, symo, i, j, trig_subs=True, sep_const=False):
    k = robo.common_root(i, j)
    chain1 = robo.chain(j, k)
    chain2 = robo.chain(i, k)
    chain2.reverse()
    complete_chain = (chain1 + chain2 + [None])
    T_out = {(j, j): eye(4)}
    T_res = eye(4)
    T = eye(4)
    for indx, x in enumerate(complete_chain[:-1]):
        inverted = indx >= len(chain1)
        T = _transform(robo, x, inverted) * T
        if trig_subs:
            for ang, name in robo.get_angles(x):
                symo.trig_replace(T, ang, name)
        T = T.expand()
        T = T.applyfunc(symo.CS12_simp)
        x_next = complete_chain[indx + 1]
        if inverted:
            t_name = (x, j)
        else:
            t_name = (robo.ant[x], j)
        T_out[t_name] = T * T_res
        if robo.paral(x, x_next):
            continue
        T_res = T * T_res
        T = eye(4)
    return T_out


def _dgm_right(robo, symo, i, j, trig_subs=True, sep_const=False):
    k = robo.common_root(i, j)
    chain1 = robo.chain(i, k)
    chain2 = robo.chain(j, k)
    chain2.reverse()
    complete_chain = (chain1 + chain2 + [None])
    T_out = {(i, i): eye(4)}
    T_res = eye(4)
    T = eye(4)
    for indx, x in enumerate(complete_chain[:-1]):
        inverted = indx >= len(chain1)
        T = T * _transform(robo, x, inverted)
        if trig_subs:
            for ang, name in robo.get_angles(x):
                symo.trig_replace(T, ang, name)
        T = T.expand()
        T = T.applyfunc(symo.CS12_simp)
        x_next = complete_chain[indx + 1]
        if inverted:
            t_name = (i, x)
        else:
            t_name = (i, robo.ant[x])
        T_out[t_name] = T_res * T
        if robo.paral(x, x_next):
            continue
        T_res = T_res * T
        T = eye(4)
    return T_out


def _dgm_one(robo, symo, i, j, fast_form=True,
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
        T = _transform(robo, x, inverted) * T
        if trig_subs:
            for ang, name in robo.get_angles(x):
                symo.trig_replace(T, ang, name)
        if fast_form:
            T = T.expand()
        T = T.applyfunc(symo.CS12_simp)
        x_next = complete_chain[indx + 1]
        if robo.paral(x, x_next):    # false if x_next is None
            continue
        T_res = T * T_res
        T = eye(4)
        if fast_form:
            _dgm_rename(robo, symo, T_res, x, i, j, inverted, forced)
    if not fast_form and forced:
        _dgm_rename(robo, symo, T_res, x, i, j, inverted, forced)
    return T_res


def _dgm_rename(robo, symo, T_res, x, i, j, inverted, forced):
    if inverted:
        name = trans_name(robo, x, j)
        forced_now = x == i
    else:
        name = trans_name(robo, robo.ant[x], j)
        forced_now = robo.ant[x] == i
    symo.mat_replace(T_res, name, forced=forced and forced_now, skip=1)


#TODO: implemet returning all the matrices
#    sep_const: bool, optional
#        If True, transform will be represented as a tuple of
#        form (Cpref, T, Cpost) where the requirad transform is
#        represented by Cpref*T*Cpost and Cpref and Cpost contain no
#        joint variables, just constant values. Only for 'left' and 'right'
def dgm(robo, symo, i, j, key='one', fast_form=True, forced=False,
        trig_subs=True):
    """must be the final DGM function

     Parameters
    ==========
    symo: Symoro
        Instance of Symoro. All the substitutions will
        be put into symo.sydi
    i: int
        To-frame index.
    j: int
        From-frame index.
    key: {'one','left','right'}
        Defines whether return just one transform or all the chain
        with multiplication from left and right
    fast_form: bool, optional
        If False, result will be in unfolded mode (triginimetric
        substitutions only)
    forced: bool, optional
        If True, all the symbols of the last transformation
        matrix will be rplaced, aplicable only if fast_form is True
        for key='one'
    trig_subs: bool, optional
        If True, all the sin(x) and cos(x) will be replaced by symbols
        SX and CX with adding them to the dictionary

    """
    if key == 'left':
        return _dgm_left(robo, symo, i, j, trig_subs)
    elif key == 'right':
        return _dgm_right(robo, symo, i, j, trig_subs)
    else:
        return _dgm_one(robo, symo, i, j, fast_form, forced, trig_subs)


def _rot(axis='z', th=0):
    """Rotation matrix about axis

    Parameters
    ==========
    axis: {'x', 'y', 'z'}
        Rotation axis
    th: var
        Rotation angle

    Returns
    =======
    rot: Matrix 3x3
    """
    if axis == 'x':
        return Matrix([[1, 0, 0],
                       [0, cos(th), - sin(th)],
                       [0, sin(th), cos(th)]])
    elif axis == 'y':
        return Matrix([[cos(th), 0, sin(th)],
                       [0, 1, 0],
                       [-sin(th), 0, cos(th)]])
    else:
        return Matrix([[cos(th), - sin(th), 0],
                       [sin(th), cos(th), 0],
                       [0, 0, 1]])


def _trans_vect(axis='z', p=0):
    """Translation vector along axis

    Parameters
    ==========
    axis: {'x', 'y', 'z'}
        Translation axis
    p: var
        Translation distance

    Returns
    =======
    v: Matrix 3x1
    """
    axis_dict = {'x': 0, 'y': 1, 'z': 2}
    v = zeros(3, 1)
    v[axis_dict[axis]] = p
    return v


def _trans(axis='z', p=0):
    """Translation matrix along axis

    Parameters
    ==========
    axis: {'x', 'y', 'z'}
        Translation axis
    p: var
        Translation distance

    Returns
    =======
    trans: Matrix 4x4
    """
    return Matrix([eye(3).row_join(_trans_vect(axis, p)),
                   [0, 0, 0, 1]])


def _rot_trans(axis='z', th=0, p=0):
    """Transformation matrix with rotation about and
    translation along axis

    Parameters
    ==========
    axis: {'x', 'y', 'z'}
        Transformation axis
    p: var
        Translation distance
    th: var
        Rotation angle

    Returns
    =======
    rot_trans: Matrix 4x4
    """
    return Matrix([_rot(axis, th).row_join(_trans_vect(axis, p)),
                   [0, 0, 0, 1]])


def compute_rot_trans(robo, symo):
    #init transformation
    antRj = Init.init_mat(robo)
    antPj = Init.init_vec(robo)
    for j in xrange(robo.NL):
        compute_transform(robo, symo, j, antRj, antPj)
    return antRj, antPj


#TODO: validate for different structures
def direct_geometric_fast(robo, i, j):
    """Computes trensformation matrix iTj.

    Parameters
    ==========
    robo: Robot
        Instance of robot description container
    i: int
        the to-frame
    j: int
        the from-frame

    Returns
    =======
    symo: Symoro
        Instance that contains all the relations of the computed model
    """
    symo = Symoro()
    symo.file_open(robo, 'fgm')
    symo.write_params_table(robo, 'Direct Geometrix model')
    dgm(robo, symo, i, j, fast_form=True, forced=True)
    symo.file_close()
    return symo


def direct_geometric(robo, frames, trig_subs):
    """Computes trensformation matrix iTj.

    Parameters
    ==========
    robo: Robot
        Instance of robot description container
    frames: list of tuples of type (i,j)
        Defines list of required transformation matrices iTj
    trig_subs: bool, optional
        If True, all the sin(x) and cos(x) will be replaced by symbols
        SX and CX with adding them to the dictionary

    Returns
    =======
    symo: Symoro
        Instance that contains all the relations of the computed model
    """
    symo = Symoro()
    symo.file_open(robo, 'trm')
    symo.write_params_table(robo, 'Direct Geometrix model')
    for i, j in frames:
        symo.write_line('Tramsformation matrix %s T %s' % (i, j))
        print dgm(robo, symo, i, j, fast_form=False,
                  forced=True, trig_subs=trig_subs)
        symo.write_line()
    symo.file_close()
    return symo
