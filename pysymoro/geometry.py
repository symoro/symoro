#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module of SYMORO package computes the geometric models.
"""

from sympy import Matrix, zeros, eye, sin, cos
from copy import copy

from symoroutils import symbolmgr
from symoroutils import tools
from symoroutils.paramsinit import ParamsInit


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


class CompTransf:
    def __init__(self, transform_type, axis, val, i=0, j=0, name=''):
        self.type = transform_type
        self.axis = axis
        self.val = val
        self.i = i
        self.j = j
        self.name = "%s%s" % (name, max(i, j))

    def __str__(self):
        axes = ('x', 'y', 'z')
        trans_type = ('rot', 'trans')
        return "%s(%s, %s)" % (trans_type[self.type],
                               axes[self.axis], self.val)

    def __repr__(self):
        return str(self)

    def matrix(self):
        """
        Homogeneous transformation matrix
        """
        if self.type == 0:
            return _rot_trans(self.axis, th=self.val)
        elif self.type == 1:
            return _rot_trans(self.axis, p=self.val)

    def matrix_inv(self):
        """
        Homogeneous transformation matrix (inverted)
        """
        if self.type == 0:
            return _rot_trans(self.axis, th=-self.val)
        elif self.type == 1:
            return _rot_trans(self.axis, p=-self.val)

    def rot(self):
        if self.type == 0:
            return _rot(self.axis, self.val)
        elif self.type == 1:
            return eye(3)

    def trans(self):
        if self.type == 1:
            return _trans_vec(self.axis, self.val)
        else:
            return zeros(3, 1)


class TransConvolve:
    def __init__(self, symo=None, trig_subs=False, simplify=True):
        self.rot = CompTransf(0, 0, 0)
        self.rot_mat = eye(3)
        self.trans = zeros(3, 1)
        self.symo = symo
        self.trig_subs = trig_subs and symo is not None
        self.T_tmp = eye(4)
        self.simplify = simplify

    def process(self, tr):
        if tr.type == 0:  # rotation
            if self.rot.axis == tr.axis and self.simplify:
                self.rot.val += tr.val
                self.rot.name += tr.name
            else:  # translation
                self.rot_mat = self.rot_mat * self.rot.rot()
                if self.trig_subs:
                    self.symo.trig_replace(self.rot_mat, self.rot.val,
                                           self.rot.name)
                self.rot = copy(tr)
        elif tr.type == 1:
            self.trans += self.rot_mat * self.rot.rot() * tr.trans()
            if self.trig_subs:
                self.symo.trig_replace(self.trans, self.rot.val,
                                       self.rot.name)

    def process_left(self, tr):
        if tr.type == 0:  # rotation
            self.trans = tr.rot() * self.trans
            if self.trig_subs:
                self.symo.trig_replace(self.trans, tr.val, tr.name)
            if self.rot.axis == tr.axis:
                self.rot.val += tr.val
                self.rot.name += tr.name
            else:  # translation
                self.rot_mat = self.rot.rot() * self.rot_mat
                if self.trig_subs:
                    self.symo.trig_replace(self.rot_mat, self.rot.val,
                                           self.rot.name)
                self.rot = copy(tr)
        elif tr.type == 1:
            self.trans += tr.trans()

    def result(self, direction='right'):
        if direction == 'right':
            r = self.rot_mat * self.rot.rot()
        elif direction == 'left':
            r = self.rot.rot() * self.rot_mat
        if self.trig_subs:
            self.symo.trig_replace(self.trans, self.rot.val, self.rot.name)
            self.symo.trig_replace(r, self.rot.val, self.rot.name)
        return Matrix([r.row_join(self.trans),
                      [0, 0, 0, 1]])


def transform_list(robo, i, j):
    """
    Computes the chain of transformations for iTj
    """
    _x = 0
    _z = 2
    k = robo.common_root(i, j)
    chain1 = robo.chain(i, k)
    chain2 = robo.chain(j, k)
    chain2.reverse()
    tr_list = []
    for indx in chain1:
        ant = robo.ant[indx]
        tr_list.append(CompTransf(0, _z, -robo.theta[indx], indx, ant))
        tr_list.append(CompTransf(1, _z, -robo.r[indx], indx, ant))
        tr_list.append(CompTransf(0, _x, -robo.alpha[indx], indx, ant, 'A'))
        tr_list.append(CompTransf(1, _x, -robo.d[indx], indx, ant))
        tr_list.append(CompTransf(0, _z, -robo.gamma[indx], indx, ant, 'G'))
        tr_list.append(CompTransf(1, _z, -robo.b[indx], indx, ant))
    for indx in chain2:
        ant = robo.ant[indx]
        tr_list.append(CompTransf(0, _z, robo.gamma[indx], ant, indx, 'G'))
        tr_list.append(CompTransf(1, _z, robo.b[indx], ant, indx))
        tr_list.append(CompTransf(0, _x, robo.alpha[indx], ant, indx, 'A'))
        tr_list.append(CompTransf(1, _x, robo.d[indx], ant, indx))
        tr_list.append(CompTransf(0, _z, robo.theta[indx], ant, indx))
        tr_list.append(CompTransf(1, _z, robo.r[indx],  ant, indx))
    return [tr for tr in tr_list if tr.val != 0]


def to_matrix_fast(symo, tr_list, forced=False):
    conv = TransConvolve(symo, trig_subs=True)
    T = eye(4)
    i = tr_list[0].i
    j = tr_list[0].j
    for tr in tr_list:
        if tr.j != j:
            T *= conv.result()
            symo.mat_replace(T, 'T%sT%s' % (i, j), skip=1)
            j = tr.j
            conv = TransConvolve(symo, trig_subs=True)
        conv.process(tr)
    T *= conv.result()
    symo.mat_replace(T, 'T%sT%s' % (i, j), skip=1, forced=forced)
    return T


def to_matrix(tr_list, symo=None, trig_subs=False, simplify=True):
    conv = TransConvolve(symo, trig_subs, simplify)
    for tr in tr_list:
        conv.process(tr)
    return conv.result()


def to_matrices_right(tr_list, symo=None, trig_subs=False):
    conv = TransConvolve(symo, trig_subs)
    j = tr_list[0].j
    i = tr_list[0].i
    res = {(i, i): eye(4)}
    for tr in tr_list:
        if tr.j != j:
            res[i, j] = conv.result()
            j = tr.j
        conv.process(tr)
    res[i, j] = conv.result()
    return res


def to_matrices_left(tr_list, symo=None, trig_subs=False):
    conv = TransConvolve(symo, trig_subs)
    j = tr_list[-1].j
    i = tr_list[-1].i
    res = {(j, j): eye(4)}
    for tr in reversed(tr_list):
        if tr.i != i:
            res[i, j] = conv.result('left')
            i = tr.i
        conv.process_left(tr)
    res[i, j] = conv.result('left')
    return res


def dgm(robo, symo, i, j, key='one', fast_form=True,
        trig_subs=True, forced=False):
    """must be the final DGM function

    Parameters
    ==========
    symo: symbolmgr.SymbolManager
        Instance of symbolmgr.SymbolManager. All the substitutions will
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
    trig_subs: bool, optional
        If True, all the sin(x) and cos(x) will be replaced by symbols
        SX and CX with adding them to the dictionary

    """
    if i == j:
        if key == 'one':
            return eye(4)
        else:
            return {(i, i): eye(4)}
    tr_list = transform_list(robo, i, j)
    if key == 'one' and fast_form:
        return to_matrix_fast(symo, tr_list, forced)
    else:
        if key == 'left':
            return to_matrices_left(tr_list, symo, trig_subs)
        elif key == 'right':
            return to_matrices_right(tr_list, symo, trig_subs)
        elif key == 'one':
            return to_matrix(tr_list, symo, trig_subs)


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
        R1 = _rot_trans(2, robo.gamma[j], robo.b[j])
        R2 = _rot_trans(0, robo.alpha[j], robo.d[j])
        R3 = _rot_trans(2, robo.theta[j], robo.r[j])
        return R1*R2*R3
    else:
        R1 = _rot_trans(2, -robo.gamma[j], -robo.b[j])
        R2 = _rot_trans(0, -robo.alpha[j], -robo.d[j])
        R3 = _rot_trans(2, -robo.theta[j], -robo.r[j])
        return R3*R2*R1


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
    ET = symo.mat_replace(-jRant*tools.skew(antPj[j]), 'JPR', j)
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3, 3).row_join(jRant)]))


def _rot(axis=2, th=0):
    """Rotation matrix about axis

    Parameters
    ==========
    axis: {0, 1, 2}
        Rotation axis
    th: var
        Rotation angle

    Returns
    =======
    rot: Matrix 3x3
    """
    assert axis in {0, 1, 2}
    if axis == 0:
        return Matrix([[1, 0, 0],
                       [0, cos(th), - sin(th)],
                       [0, sin(th), cos(th)]])
    elif axis == 1:
        return Matrix([[cos(th), 0, sin(th)],
                       [0, 1, 0],
                       [-sin(th), 0, cos(th)]])
    else:
        return Matrix([[cos(th), - sin(th), 0],
                       [sin(th), cos(th), 0],
                       [0, 0, 1]])


def _trans_vec(axis=2, p=0):
    """Translation vector along axis

    Parameters
    ==========
    axis: {0, 1, 2}
        Translation axis
    p: var
        Translation distance

    Returns
    =======
    v: Matrix 3x1
    """
    assert axis in {0, 1, 2}
    v = zeros(3, 1)
    v[axis] = p
    return v


def _rot_trans(axis=2, th=0, p=0):
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
    assert axis in {0, 1, 2}
    return Matrix([_rot(axis, th).row_join(_trans_vec(axis, p)),
                   [0, 0, 0, 1]])


def compute_rot_trans(robo, symo):
    #init transformation
    antRj = ParamsInit.init_mat(robo)
    antPj = ParamsInit.init_vec(robo)
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
    symo: symbolmgr.SymbolManager
        Instance that contains all the relations of the computed model
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'fgm')
    symo.write_params_table(robo, 'Direct Geometric model')
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
    symo: symbolmgr.SymbolManager
        Instance that contains all the relations of the computed model
    """
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'trm')
    symo.write_params_table(robo, 'Direct Geometric model')
    for i, j in frames:
        symo.write_line('Tramsformation matrix %s T %s' % (i, j))
        T = dgm(robo, symo, i, j, fast_form=False, trig_subs=trig_subs)
        symo.mat_replace(T, 'T%sT%s' % (i, j), forced=True, skip=1)
        symo.write_line()
    symo.file_close()
    return symo
