# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package computes the geometric models.
"""


from sympy import Matrix, zeros, eye, sin, cos


from symoroutils.tools import ZERO, skew, simplify
from symoroutils.paramsinit import ParamsInit


Z_AXIS = Matrix([0, 0, 1])


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


class Transform(Matrix):
    def __init__(self, *args):
        Matrix.__init__(*args)
        assert self.shape == (4, 4)

    @classmethod
    def create(self, axis=2, th=0, p=0):
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
        assert axis in (0, 1, 2)
        R = _rot(axis, th)
        P = _trans_vec(axis, p)
        return self.from_RP(R, P)

    @classmethod
    def from_RP(self, R, P):
        assert R.shape == (3, 3)
        assert P.shape == (3, 1)
        return Transform([R.row_join(P), [0, 0, 0, 1]])

    @classmethod
    def frame(self, robo, j):
        """Transform matrix between frames j and ant[j]

        Parameters
        ==========
        j: int
            Frame index.

        Returns
        =======
        transform: Matrix 4x4
            Transformation matrix. If invert is True then j_T_ant,
            else ant_T_j.
        """
        T1 = Transform.create(2, robo.gamma[j], robo.b[j])
        T2 = Transform.create(0, robo.alpha[j], robo.d[j])
        T3 = Transform.create(2, robo.theta[j], robo.r[j])
        return T1*T2*T3

    @classmethod
    def frame_inv(self, robo, j):
        """Transform matrix between frames j and ant[j]

        Parameters
        ==========
        j: int
            Frame index.

        Returns
        =======
        transform: Matrix 4x4
            Transformation matrix. If invert is True then j_T_ant,
            else ant_T_j.
        """
        T1 = Transform.create(2, -robo.gamma[j], -robo.b[j])
        T2 = Transform.create(0, -robo.alpha[j], -robo.d[j])
        T3 = Transform.create(2, -robo.theta[j], -robo.r[j])
        return T3*T2*T1

    @property
    def sna(self):
        """Extracts the s, n, a vector basis of rotation 3x3 matrix
        from 4x4 transformation matrix

        Returns
        =======
        s: Matrix 3x1
        n: Matrix 3x1
        a: Matrix 3x1
        """
        R = Matrix(self.R)
        return R.col(0), R.col(1), R.col(2)

    @property
    def R(self):
        """Extracts rotation 3x3 matrix from 4x4 transformation matrix

        Returns
        =======
        get_r: Matrix 3x3
        """
        return Matrix(self[:3, :3])

    @property
    def P(self):
        """Extracts translation vector from 4x4 transformation matrix

        Returns
        =======
        get_p: Matrix 3x1
        """
        return Matrix(self[:3, 3])

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
    def kPj(self, robo, antPj, antRj, k, chainj):
        T = eye(4)
        for i in chainj:
            if i > k:
                kTj = eye(4)
                kTj[0:3, 0:3] = antRj[i]
                kTj[0:3, 3:] = antPj[i]
                T = kTj * T
        Px = T[0, 3]
        Py = T[1, 3]
        Pz = T[2, 3]
        return Px, Py, Pz

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

    @property
    def along_z(self):
        return self.R.col(2) == Z_AXIS


class DGM:
    @classmethod
    def compute(self, robo, symo, i, j, fast_form=False, trig_subs=False):
        """Function that computes the transformation matrix between
        frames i and j. That is iTj.

        Parameters
        ==========
        symo: symbolmgr.SymbolManager
            Instance of symbolmgr.SymbolManager. All the substitutions will
            be put into symo.sydi
        i: int
            To-frame index.
        j: int
            From-frame index.
        fast_form: bool, optional
            If False, result will be in unfolded mode (triginimetric
            substitutions only)
        trig_subs: bool, optional
            If True, all the sin(x) and cos(x) will be replaced by symbols
            SX and CX with adding them to the dictionary
        forced; bool, optional
            If True, the final matrix elements will be replaced by symbols
            like TiTjab even if the are already simple or even constants.
        """
        dgm_acc = DGM(robo, symo, fast_form=fast_form, trig_subs=trig_subs)
        return dgm_acc.transform(i, j)

    @classmethod
    def compute_left(self, robo, symo, i, j, trig_subs=False):
        """This function computes all the transformation matrices
        kTj for k in [i..j]

        Parameters
        ==========
        symo: symbolmgr.SymbolManager
            Instance of symbolmgr.SymbolManager. All the substitutions will
            be put into symo.sydi
        i: int
            To-frame index.
        j: int
            From-frame index.
        trig_subs: bool, optional
            If True, all the sin(x) and cos(x) will be replaced by symbols
            SX and CX with adding them to the dictionary
        """
        dgm_acc = DGM(robo, symo, trig_subs=trig_subs)
        res = dict()
        for k in robo.loop_chain(i, j):
            res[k, j] = dgm_acc.transform(k, j)
        return res

    @classmethod
    def compute_right(self, robo, symo, i, j, trig_subs=False):
        """This function computes all the transformation matrices
        iTk for k in [i..j]

        Parameters
        ==========
        symo: symbolmgr.SymbolManager
            Instance of symbolmgr.SymbolManager. All the substitutions will
            be put into symo.sydi
        i: int
            To-frame index.
        j: int
            From-frame index.
        trig_subs: bool, optional
            If True, all the sin(x) and cos(x) will be replaced by symbols
            SX and CX with adding them to the dictionary
        """
        dgm_acc = DGM(robo, symo, trig_subs=trig_subs)
        res = dict()
        for k in robo.loop_chain(i, j):
            res[i, k] = dgm_acc.transform(i, k)
        return res

    def __init__(self, robo, symo=None, trig_subs=False, fast_form=False):
        self.symo = symo
        self.robo = robo
        self.cash = dict()
        self.inv_cash = dict()
        if symo is not None:
            self.trig_subs = trig_subs or fast_form
            self.fast_form = fast_form

    def block_chain(self, chain):
        """Computes blocks of parallel axes
        For example output for RX90 is [6, 5, 4, [3, 2], 1]
        """
        res = []
        paral_block = []
        for k in chain:
            paral_block.append(k)
            if self.robo.alpha[k] != ZERO or k == chain[-1]:
                if len(paral_block) > 1:
                    res.append(paral_block)
                else:
                    res.append(k)
                paral_block = []
        return res

    def multitransform(self, block):
        if isinstance(block, list):
            j = block[0]
            key = (j, block[-1])
            if key in self.cash:
                return self.cash[key]
            else:
                T = Transform.frame(self.robo, block[0])
                for i in block[1:]:
                    T = simplify(Transform.frame(self.robo, i) * T)
                    self.cash[(j, i)] = T
        else:
            if block in self.cash:
                return self.cash[block]
            else:
                T = Transform.frame(self.robo, block)
                self.cash[block] = T
        return T

    def multitransform_inv(self, block):
        if isinstance(block, list):
            j = block[0]
            key = (j, block[-1])
            if key in self.inv_cash:
                return self.inv_cash[key]
            else:
                T = Transform.frame_inv(self.robo, j)
                for i in block[1:]:
                    T = simplify(T * Transform.frame_inv(self.robo, i))
                    self.inv_cash[(j, i)] = T
        else:
            if block in self.inv_cash:
                return self.inv_cash[block]
            else:
                T = Transform.frame_inv(self.robo, block)
                self.inv_cash[block] = T
        return T

    def transform(self, i, j):
        if i == j:
            return Transform(eye(4))
        k = self.robo.common_root(i, j)
        chain1 = self.robo.chain(i, k)
        chain2 = self.robo.chain(j, k)
        b_chain1 = self.block_chain(chain1)
        b_chain2 = self.block_chain(chain2)
        b_chain2.reverse()
        T_res = Transform(eye(4))
        for idx, k in enumerate(b_chain1):
            T_res = T_res*self.multitransform_inv(k)
            if self.fast_form:
                self.symo.trig_replace(T_res)
                self.symo.mat_replace(T_res, 'TI%s' % idx)
        for idx, k in enumerate(b_chain2):
            T_res = T_res*self.multitransform(k)
            if self.fast_form:
                self.symo.trig_replace(T_res)
                self.symo.mat_replace(T_res, 'TD%s' % idx)
        if self.trig_subs and not self.fast_form:
            self.symo.trig_replace(T_res)
        return T_res


def compute_transform(robo, symo, j, antRj, antPj):
    """Internal function. Computes rotation matrix and translation vector
    of ant_T_j homogenuous transform. Does the trigonometric subsctitution
    and saves the symbols into symo.sydi

    Notes
    =====
    antPj and antRj are the output parameters
    """
    antTj = Transform.frame(robo, j)
    for angle, name in robo.get_angles(j):
        antTj = symo.trig_replace(antTj)
    antRj[j] = symo.mat_replace(antTj.R, 'A', j)
    antPj[j] = symo.mat_replace(antTj.P, 'L', j)


def compute_screw_transform(robo, symo, j, antRj, antPj, jTant):
    """Internal function. Computes the screw transformation matrix
    between ant[j] and j frames.

    Notes
    =====
    jTant is an output parameter
    """
    jRant = antRj[j].T
    ET = symo.mat_replace(-jRant*skew(antPj[j]), 'JPR', j)
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3, 3).row_join(jRant)]))


def compute_rot_trans(robo, symo):
    #init transformation
    antRj = ParamsInit.init_mat(robo)
    antPj = ParamsInit.init_vec(robo)
    for j in xrange(robo.NL):
        compute_transform(robo, symo, j, antRj, antPj)
    return antRj, antPj


##TODO: bring to the interface file
#def direct_geometric_fast(robo, i, j):
#    """Computes trensformation matrix iTj.
#
#    Parameters
#    ==========
#    robo: Robot
#        Instance of robot description container
#    i: int
#        the to-frame
#    j: int
#        the from-frame
#
#    Returns
#    =======
#    symo: symbolmgr.SymbolManager
#        Instance that contains all the relations of the computed model
#    """
#    symo = symbolmgr.SymbolManager()
#    symo.file_open(robo, 'fgm')
#    symo.write_params_table(robo, 'Direct Geometric model')
#    T = DGM.compute(robo, symo, i, j, fast_form=True)
#    symo.mat_replace(T, 'T%sT%s' % (i, j), forced=True, skip=1)
#    symo.file_close()
#    return symo
#
#
#def direct_geometric(robo, frames, trig_subs):
#    """Computes trensformation matrix iTj.
#
#    Parameters
#    ==========
#    robo: Robot
#        Instance of robot description container
#    frames: list of tuples of type (i,j)
#        Defines list of required transformation matrices iTj
#    trig_subs: bool, optional
#        If True, all the sin(x) and cos(x) will be replaced by symbols
#        SX and CX with adding them to the dictionary
#
#    Returns
#    =======
#    symo: symbolmgr.SymbolManager
#        Instance that contains all the relations of the computed model
#    """
#    symo = symbolmgr.SymbolManager()
#    symo.file_open(robo, 'trm')
#    symo.write_params_table(robo, 'Direct Geometric model')
#    dgm = DGM(robo, symo, trig_subs=trig_subs)
#    for i, j in frames:
#        symo.write_line('Tramsformation matrix %s T %s' % (i, j))
#        T = dgm.transform(i, j)
#        symo.mat_replace(T, 'T%sT%s' % (i, j), forced=True, skip=1)
#        symo.write_line()
#    symo.file_close()
#    return symo


