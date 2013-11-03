#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This module of SYMORO package provides kinematic models' computation.

The core symbolic library is sympy.

Needed modules : symoro.py, geometry.py

ECN - ARIA1 2013
"""


from sympy import Matrix, zeros

try:
    from pysymoro.core.symoro import Symoro, Init, hat
    from pysymoro.core.symoro import FAIL, ZERO
    from pysymoro.core.geometry import dgm, Transform, compute_rot_trans
except ImportError:
    from core.symoro import Symoro, Init, hat
    from core.symoro import FAIL, ZERO
    from core.geometry import dgm, Transform, compute_rot_trans


TERMINAL = 0
ROOT = 1


def compute_omega(robo, symo, j, antRj, w, wi):
    """Internal function. Computes angular velocity of jth frame and
    projection of the antecedent frame's angular velocity

    Notes
    =====
    w, wi, U, vdot are the output parameters
    """
    jRant = antRj[j].T
    qdj = Matrix([0, 0, robo.qdot[j]])
    wi[j], w[j] = omega_ij(robo, symo, j, jRant, w, qdj)


def omega_ij(robo, symo, j, jRant, w, qdj, forced=False):
    wi = symo.mat_replace(jRant*w[robo.ant[j]], 'WI', j)
    wj = symo.mat_replace(wi + (1 - robo.sigma[j])*qdj, 'W', j, forced)
    return wi, wj


def compute_twist(robo, symo, j, antRj, antPj, w, wdot, U, vdot, forced=False):
    """Internal function. Computes angular velocity, auxiliary U matrix and
    linear and angular accelerations.

    Notes
    =====
    w, wdot, U, vdot are the output parameters
    """
    jRant = antRj[j].T
    qdj = Matrix([0, 0, robo.qdot[j]])
    qddj = Matrix([0, 0, robo.qddot[j]])
    wi, w[j] = omega_ij(robo, symo, j, jRant, w, qdj, forced)
    self_rotation = (1 - robo.sigma[j])*(qddj + hat(wi)*qdj)
    wdot[j] = jRant*wdot[robo.ant[j]] + self_rotation
    symo.mat_replace(wdot[j], 'WP', j, forced)
    DV = Init.product_combinations(w[j])
    symo.mat_replace(DV, 'DV', j)
    hatw_hatw = Matrix([[-DV[3]-DV[5], DV[1], DV[2]],
                        [DV[1], -DV[5]-DV[0], DV[4]],
                        [DV[2], DV[4], -DV[3]-DV[0]]])
    U[j] = hatw_hatw + hat(wdot[j])
    symo.mat_replace(U[j], 'U', j)
    vsp = vdot[robo.ant[j]] + U[robo.ant[j]]*antPj[j]
    symo.mat_replace(vsp, 'VSP', j)
    vdot[j] = jRant*vsp + robo.sigma[j]*(qddj + 2*hat(wi)*qdj)
    symo.mat_replace(vdot[j], 'VP', j, forced)


def _jac(robo, symo, n, i, j, chain=None, forced=False, trig_subs=False):
    """
    Computes jacobian of frame n (with origin On in Oj) projected to frame i
    """
#    symo.write_geom_param(robo, 'Jacobian')
    M = []
    if chain is None:
        chain = reversed(robo.chain(n))
    for k in chain:
        kTj = dgm(robo, symo, k, j, fast_form=False, trig_subs=trig_subs)
        iTk = dgm(robo, symo, i, k, fast_form=False, trig_subs=trig_subs)
        isk, ink, iak = Transform.sna(iTk)
        sigm = robo.sigma[k]
        if sigm == 1:
            dvdq = iak
            J_col = dvdq.col_join(Matrix([0, 0, 0]))
        elif sigm == 0:
            dvdq = kTj[0, 3]*ink-kTj[1, 3]*isk
            J_col = dvdq.col_join(iak)
        else:
            J_col = Matrix([0, 0, 0, 0, 0, 0])
        M.append(J_col.T)
    Jac = Matrix(M).T
    Jac = Jac.applyfunc(symo.simp)
    iTj = dgm(robo, symo, i, j, fast_form=False, trig_subs=trig_subs)
    jTn = dgm(robo, symo, j, n, fast_form=False, trig_subs=trig_subs)
    iRj = Transform.R(iTj)
    jPn = Transform.P(jTn)
    L = -hat(iRj*jPn)
    if forced:
#        symo.write
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
        J, L = _jac(robo, symo, n, i, j, False)
    if not J.is_square:
        J = _make_square(J)
    det = J.det()
    det = symo.simp(det)
    return det


def extend_W(J, r, W, indx, chain):
    row = []
    for e in indx:
        if e in chain:
            row.append(J[r, chain.index(e)])
        else:
            row.append(0)
    W.append(row)


def _kinematic_loop_constraints(robo, symo, proj=None):
    B = robo.NJ - robo.NL
    if B == 0:
        print 'There are no loops'
        return FAIL
    indx_c = range(robo.NL, robo.NL + B)
    indx_a, indx_p = [], []
    for i in xrange(1, robo.NL):
        if robo.mu[i] == 1:
            indx_a.append(i)
        else:
            indx_p.append(i)
    W_a, W_p, W_ac, W_pc, W_c = [], [], [], [], []
    for indx, (i, j) in enumerate(robo.loop_terminals):
        # i - cut joint, j - fixed joint
        k = robo.common_root(i, j)
        chi = robo.chain(i, k)
        chj = robo.chain(j, k)
        if proj is not None and len(proj) == B and proj[indx] == TERMINAL:
            Ji, L = _jac(robo, symo, i, i, i, chi)
            Jj, L = _jac(robo, symo, j, j, j, chj)
        else:
            Ji, L = _jac(robo, symo, k, i, i, chi)
            Jj, L = _jac(robo, symo, k, j, j, chj)
        chi.extend(chj)
        J = Ji.row_join(-Jj)
        for row in xrange(6):
            if all(J[row, col] == ZERO for col in xrange(len(chi))):
                continue
            elif J[row, chi.index(i)] == ZERO:
                extend_W(J, row, W_a, indx_a, chi)
                extend_W(J, row, W_p, indx_p, chi)
            else:
                extend_W(J, row, W_ac, indx_a, chi)
                extend_W(J, row, W_pc, indx_p, chi)
                extend_W(J, row, W_c, indx_c, chi)
    W_a, W_p = Matrix(W_a), Matrix(W_p)
    W_ac, W_pc, W_c = Matrix(W_ac), Matrix(W_pc), Matrix(W_c)
    print W_a, W_p, W_ac, W_pc, W_c
    return W_a, W_p, W_ac, W_pc, W_c


def compute_speeds_accelerations(robo, symo, antRj=None, antPj=None):
    """Internal function. Computes speeds and accelerations usitn

    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    symo : Symoro
        Instance of symbolic manager
    """

    #init velocities and accelerations
    w = Init.init_w(robo)
    wdot, vdot = Init.init_wv_dot(robo)
    #init transformation
    if antRj is None or antPj is None:
        antRj, antPj = compute_rot_trans(robo, symo)
        forced = True
    else:
        forced = False
    #init auxilary matrix
    U = Init.init_U(robo)

    for j in xrange(1, robo.NL):
        compute_twist(robo, symo, j, antRj, antPj, w, wdot, U, vdot, forced)
    return w, wdot, vdot, U


def speeds_accelerations(robo):
    symo = Symoro(None)
    symo.file_open(robo, 'acc')
    symo.write_params_table(robo, 'Speeds and accelerations')
    compute_speeds_accelerations(robo, symo)
    symo.file_close()
    return symo


def jacobian(robo, n, i, j):
    symo = Symoro()
    symo.file_open(robo, 'jac')
    title = "Jacobian matrix for frame %s\n"
    title += "Projection frame %s, intermediat frame %s"
    symo.write_params_table(robo, title % (n, i, j))
    _jac(robo, symo, n, i, j, forced=True)
    symo.file_close()
    return symo


def jacobian_determinant(robo, n, i, j, rows, cols):
    symo = Symoro(None)
    J, L = _jac(robo, symo, n, i, j, trig_subs=False)
    J_reduced = zeros(len(rows), len(cols))
    for i, i_old in enumerate(rows):
        for j, j_old in enumerate(cols):
            J_reduced[i, j] = J[i_old, j_old]
    symo.file_open(robo, 'det')
    symo.write_params_table(robo, 'Jacobian determinant for frame %s' % n)
    symo.write_line(_jac_det(robo, symo, J=J_reduced))
    symo.file_close()
    return symo


def kinematic_constraints(robo):
    symo = Symoro(None)
    res = _kinematic_loop_constraints(robo, symo)
    if res == FAIL:
        return FAIL
    W_a, W_p, W_ac, W_pc, W_c = res
    symo.file_open(robo, 'ckel')
    symo.write_params_table(robo, 'Constraint kinematic equations of loop')
    symo.mat_replace(W_a, 'WA', forced=True)
    symo.mat_replace(W_p, 'WP', forced=True)
    symo.mat_replace(W_ac, 'WPA', forced=True)
    symo.mat_replace(W_pc, 'WPC', forced=True)
    symo.mat_replace(W_c, 'WC', forced=True)
    symo.file_close()
    return symo


#symo = Symoro()
#from symoro import Symoro, Robot
#kinematic_constraints(Robot.SR400())
##jacobian_determinant(robo, 6, range(6), range(6))
###print _jac(robo, symo, 2, 5, 5)
###print _jac_det(robo, symo, 5)
###W = kinematic_loop_constraints(robo, symo)
###print W[0]
###print W[1]
###speeds_accelerations(robo, symo)
###print _jac_inv(Robot.RX90(), symo, 2, 5, 5)
##
#def b():
#    symo = Symoro()
#    print _jac_inv(Robot.RX90(), symo, 6, 3, 3)
####from timeit import timeit
#####print timeit(a, number=10)
#####print timeit(b, number=10)
####
#import profile
###
##profile.run('b()', sort = 'cumtime')
##profile.run('b()')
#from timeit import timeit
#print timeit(b, number=1)
