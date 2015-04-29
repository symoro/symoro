# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package computes the kinematic models.
"""

#TODO think about trigsimp
from sympy import Matrix, zeros, factor

from pysymoro.geometry import DGM
from pysymoro.geometry import Z_AXIS

from symoroutils.tools import skew, ZERO, simplify, CLOSED_LOOP
from symoroutils.paramsinit import ParamsInit


TERMINAL = 0
ROOT = 1


def _omega_ij(robo, j, jRant, w, qdj):
    wi = jRant*w[robo.ant[j]]
    w[j] = wi
    if robo.sigma[j] == 0:    # revolute joint
        w[j] += qdj
    return wi, w[j]


def _omega_dot_j(robo, j, jRant, w, wi, wdot, qdj, qddj):
    wdot[j] = jRant*wdot[robo.ant[j]]
    if robo.sigma[j] == 0:    # revolute joint
        wdot[j] += (qddj + skew(wi)*qdj)
    return wdot[j]


def _v_j(robo, j, antPj, jRant, v, w, qdj, forced=False):
    ant = robo.ant[j]
    v[j] = jRant*(skew(w[ant])*antPj[j] + v[ant])
    if robo.sigma[j] == 1:     # prismatic joint
        v[j] += qdj
    return v[j]


def _v_dot_j(robo, symo, j, jRant, antPj,
             w, wi, wdot, U, vdot, qdj, qddj):
    DV = ParamsInit.product_combinations(w[j])
    symo.mat_replace(DV, 'DV', j)
    hatw_hatw = Matrix([[-DV[3]-DV[5], DV[1], DV[2]],
                        [DV[1], -DV[5]-DV[0], DV[4]],
                        [DV[2], DV[4], -DV[3]-DV[0]]])
    U[j] = hatw_hatw + skew(wdot[j])
    symo.mat_replace(U[j], 'U', j)
    vsp = vdot[robo.ant[j]] + U[robo.ant[j]]*antPj[j]
    symo.mat_replace(vsp, 'VSP', j)
    vdot[j] = jRant*vsp
    if robo.sigma[j] == 1:    # prismatic joint
        vdot[j] += qddj + 2*skew(wi)*qdj
    return vdot[j]


#TODO: figure out why it is needed
def compute_omega(robo, symo, j, antRj, w, wi):
    """Internal function. Computes angular velocity of jth frame and
    projection of the antecedent frame's angular velocity

    Notes
    =====
    w, wi are the output parameters
    """
    jRant = antRj[j].T
    qdj = Z_AXIS * robo.qdot[j]
    wi[j], w[j] = _omega_ij(robo, j, jRant, w, qdj)
    symo.mat_replace(wi[j], 'WI', j)
    symo.mat_replace(w[j], 'W', j)


def _jac(robo, symo, n, i, j, chain=None, forced=False, trig_subs=False):
    """
    Computes jacobian of frame n (with origin On in Oj) projected to frame i
    """
    #  symo.write_geom_param(robo, 'Jacobian')
    # TODO: Check projection frames, rewrite DGM call for higher efficiency
    J_col_list = []
    if chain is None:
        chain = robo.chain(n)
        chain.reverse()
    inter_chain = robo.chain(j)
    assert j in chain or j == robo.ant[chain[-1]]
    assert i in inter_chain or i == robo.ant[inter_chain[-1]]
    dgm = DGM(robo, symo)
    iTk_dict = dict()
    for k in chain:
        iTk_dict[i, k] = dgm.transform(i, k)
    iTj = iTk_dict[i, j]
    for k in chain:
        iTk = iTk_dict[i, k]
        isk, ink, iak = iTk.sna
        sigm = robo.sigma[k]
        if sigm == 1:
            dvdq = iak
            J_col = dvdq.col_join(zeros(3, 1))
        elif sigm == 0:
            dvdq = skew(iak) * (iTj.P - iTk.P)
            J_col = dvdq.col_join(iak)
        else:
            J_col = zeros(6, 1)
        J_col_list.append(J_col.T)
    Jac = Matrix(J_col_list).T
    #Jac = Jac.applyfunc(symo.simp)
    iRj = iTk_dict[i, j].R
    jTn = DGM.compute(robo, symo, j, n, trig_subs=trig_subs)
    jPn = jTn.P
    L = -skew(iRj*jPn)
    #L = L.applyfunc(symo.simp)
    if trig_subs:
        pass
    if forced:
        symo.mat_replace(Jac, 'J', '', forced)
        L = symo.mat_replace(L, 'L', '', forced)
    return Jac, L


def _make_square(J):
    if J.shape[0] > J.shape[1]:
        return J.T*J
    else:
        return J*J.T


def _jac_inv(robo, symo, n, i, j):
    J, L = _jac(robo, symo, n, i, j)
    if not J.is_square:
        J = _make_square(J)
    det = _jac_det(robo, symo, J=J)
    Jinv = J.adjugate()
    if det == ZERO:
        print 'Matrix is singular!'
    else:
        Jinv = Jinv/det
    Jinv = Jinv.applyfunc(symo.simp)
    symo.mat_replace(Jinv, 'JI', '', False)
    return Jinv


def _jac_det(robo, symo, n=1, i=1, j=1, J=None):
    if J is None:
        J, L = _jac(robo, symo, n, i, j)
    if not J.is_square:
        J = _make_square(J)
    print J
    det = J.det()
    print det
    print
    det = simplify(det, C2S2=True)
    print det
    print
    det = factor(det)
    print det
    print
    return det


def _extend_W(J, r, W, indx, chain):
    row = []
    for e in indx:
        if e in chain:
            row.append(J[r, chain.index(e)])
        else:
            row.append(0)
    W.append(row)


def kin_loop_constr(robo, symo, proj=None):
#    if robo.NJ == robo.NL:
#        return FAIL
    nloops = robo.NJ - robo.NL
    print len(proj), nloops
    assert robo.structure == CLOSED_LOOP
    assert nloops > 0
    assert proj is None or len(proj) == nloops
    indx_c = robo.indx_cut
    indx_a = robo.indx_active
    indx_p = robo.indx_passive
    print indx_p, indx_c
    W_a, W_p, W_ac, W_pc, W_c = [], [], [], [], []
    twist_rows = {}
    twist_rows_cut = {}
    row_p = 0
    row_cut = 0
    for indx, (i, j) in enumerate(robo.loop_terminals):
        # i - cut joint, j - fixed joint
        k = robo.common_root(i, j)
        chi = robo.chain(i, k)
        chj = robo.chain(j, k)
        if proj is not None and proj[indx] == TERMINAL:
            Ji, L = _jac(robo, symo, i, i, i, chi)
            Jj, L = _jac(robo, symo, j, j, j, chj)
        else:
            Ji, L = _jac(robo, symo, i, k, i, chi)
            Jj, L = _jac(robo, symo, j, k, j, chj)
        chi.extend(chj)
        J = Ji.row_join(-Jj)
        twist_rows[(i, j)] = []
        for row in xrange(6):
            if all(J[row, col] == ZERO for col in xrange(len(chi))):
                continue
            elif J[row, chi.index(i)] == ZERO:
                _extend_W(J, row, W_a, indx_a, chi)
                _extend_W(J, row, W_p, indx_p, chi)
                twist_rows[(i, j)].append((row, row_p))
                row_p += 1
            else:
                _extend_W(J, row, W_ac, indx_a, chi)
                _extend_W(J, row, W_pc, indx_p, chi)
                _extend_W(J, row, W_c, indx_c, chi)
                twist_rows_cut[(i, j)] = (row, row_cut)
                row_cut += 1
    for k, (row, row_c) in twist_rows_cut.iteritems():
        twist_rows[k].append((row, row_c + row_p))
    W_a, W_p = Matrix(W_a), Matrix(W_p)
    W_ac, W_pc, W_c = Matrix(W_ac), Matrix(W_pc), Matrix(W_c)
    return W_a, W_p, W_ac, W_pc, W_c, twist_rows


#TODO: check that gravity=True is explicit everywhere
def compute_vel_acc(robo, symo, antRj, antPj,
                    floating=True, forced=False, gravity=False):
    """Internal function. Computes speeds and accelerations usitn

    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    symo : symbolmgr.SymbolManager
        Instance of symbolic manager
    """
    #init velocities and accelerations
    w = ParamsInit.init_w(robo)
    wdot, vdot = ParamsInit.init_wv_dot(robo, gravity)
    # decide first link
    first_link = 1
    if floating or robo.is_floating or robo.is_mobile:
        first_link = 0
    #init auxilary matrix
    U = ParamsInit.init_u(robo)
    for j in xrange(first_link, robo.NL):
        if j == 0:
            w[j] = symo.mat_replace(w[j], 'W', j)
            wdot[j] = symo.mat_replace(wdot[j], 'WP', j)
            vdot[j] = symo.mat_replace(vdot[j], 'VP', j)
            dv0 = ParamsInit.product_combinations(w[j])
            symo.mat_replace(dv0, 'DV', j)
            hatw_hatw = Matrix([
                [-dv0[3]-dv0[5], dv0[1], dv0[2]],
                [dv0[1], -dv0[5]-dv0[0], dv0[4]],
                [dv0[2], dv0[4], -dv0[3]-dv0[0]]
            ])
            U[j] = hatw_hatw + skew(wdot[j])
            symo.mat_replace(U[j], 'U', j)
        else:
            jRant = antRj[j].T
            qdj = Z_AXIS * robo.qdot[j]
            qddj = Z_AXIS * robo.qddot[j]
            wi, w[j] = _omega_ij(robo, j, jRant, w, qdj)
            symo.mat_replace(w[j], 'W', j)
            symo.mat_replace(wi, 'WI', j)
            _omega_dot_j(robo, j, jRant, w, wi, wdot, qdj, qddj)
            symo.mat_replace(wdot[j], 'WP', j, forced)
            _v_dot_j(robo, symo, j, jRant, antPj,
                     w, wi, wdot, U, vdot, qdj, qddj)
            symo.mat_replace(vdot[j], 'VP', j, forced)
    return w, wdot, vdot, U

