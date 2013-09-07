"""
This module of SYMORO package provides kinematic models' computation.

The core symbolic library is sympy.

Needed modules : symoro.py, geometry.py

ECN - ARIA1 2013
"""
from math import floor
from sympy import Matrix, zeros
from symoro import Symoro, Init, hat, ZERO
from geometry import dgm, Transform, compute_rot_trans



def compute_omega(robo, symo, j, antRj, w, wi):
    """Internal function. Computes angular velocity of jth frame and
    projection of the antecedent frame's angular velocity

    Notes
    =====
    w, wi, U, vdot are the output parameters
    """
    jRant = antRj[j].T
    qdj = Matrix([0, 0, robo.qdot[j]])
    wi[j], w[j] =  omega_ij(robo, symo, j, jRant, w, qdj)

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
    wi, w[j] =  omega_ij(robo, symo, j, jRant, w, qdj, forced)
    wdot[j] = jRant*wdot[robo.ant[j]] + (1 - robo.sigma[j])*(qddj + hat(wi)*qdj)
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

def _jac(robo, symo, i, j, n, chain=None, forced=False, trig_subs=False):
    """
    Computes jacobian of frame n (with origin On in Oj) projected to frame i
    """
#    symo.write_geom_param(robo, 'Jacobian')
    M = []
    if chain == None:
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
    iTj = dgm(robo, symo, i, j, False)
    jTn = dgm(robo, symo, j, n, False)
    iRj = Transform.R(iTj)
    jPn = Transform.P(jTn)
    L = symo.mat_replace(-hat(iRj*jPn), 'L', '', forced)
    return Jac, L

def _jac_inv(robo, symo, i, j, n):
    J, L = _jac(robo, symo, i, j, n)
    if not J.is_square:
        J = J*J.T
    det = _jac_det(robo, symo, n)
    Jinv = J.adjugate()
    if det == ZERO:
        print 'Matrix is singular!'
    else:
        Jinv = Jinv/det
    Jinv = Jinv.applyfunc(symo.simp)
    symo.mat_replace(Jinv, 'JI', '', False)
    return Jinv

def _jac_det(robo, symo, n=1, J=None):
    if J == None:
        i = int(floor(n/2))
        J, L = _jac(robo, symo, i, i, n, False)
    if not J.is_square:
        J = J*J.T
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

def _kinematic_loop_constraints(robo, symo):
    B = robo.NJ - robo.NL
    if B == 0:
        print 'There are no loops'
        return
    indx_c = range(robo.NL, robo.NJ + B)
    indx_a, indx_p = [], []
    for i in xrange(1, robo.NL):
        if robo.mu[i] == 1:
            indx_a.append(i)
        else:
            indx_p.append(i)
    W_a, W_p, W_ac, W_pc, W_c = [], [], [], [], []
    for i, j in robo.get_loop_terminals():
        k = robo.common_root(i, j)
        chi = robo.chain(i, k)
        chj = robo.chain(j, k)
        Ji, L = _jac(robo, symo, k, i, i, chi)
        Jj, L = _jac(robo, symo, k, j, j, chj)
        chi.extend(chj)
        J = Ji.row_join(-Jj)
        for row in xrange(6):
            if all(J[row, col] == ZERO for col in xrange(len(chi))):
                continue
            elif J[row, chi.index(j)] == ZERO:
                extend_W(J, row, W_a, indx_a, chi)
                extend_W(J, row, W_p, indx_p, chi)
            else:
                extend_W(J, row, W_ac, indx_a, chi)
                extend_W(J, row, W_pc, indx_p, chi)
                extend_W(J, row, W_c, indx_c, chi)
    W_a, W_p = Matrix(W_a), Matrix(W_p)
    W_ac, W_pc, W_c = Matrix(W_ac), Matrix(W_pc), Matrix(W_c)
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
    if antRj == None or antPj == None:
        antRj, antPj = compute_rot_trans(robo, symo)
        forced = True
    else:
        forced = False
    #init auxilary matrix
    U = Init.init_U(robo)

    for j in xrange(1, robo.NL):
        compute_twist(robo, symo, j, antRj, antPj, w, wdot, U, vdot, forced)
    return w, wdot, vdot, U

def jacobian(robo, n, i, j):
    symo = Symoro()
    symo.file_open(robo, 'jac')
    symo.write_geom_param(robo, 'Jacobian matrix for frame %s' % n)
    symo.write_line('Projection frame %s, intermediat frame %s' % (i, j))
    symo.write_line()
    _jac(robo, symo, i, j, n)
    symo.file_out.close()
    return symo

def jacobian_determinant(robo, n, rows, cols):
    symo = Symoro(None)
    i = int(floor(n/2))
    J, L = _jac(robo, symo, i, i, n, trig_subs=False)
    J_reduced = zeros(len(rows), len(cols))
    for i, i_old in enumerate(rows):
        for j, j_old in enumerate(cols):
            J_reduced[i, j] = J[i_old, j_old]
    symo.file_open(robo, 'det')
    symo.write_geom_param(robo, 'Jacobian determinant for frame %s' % n)
    symo.write_line(_jac_det(robo, symo, J=J_reduced))
    symo.file_out.close()
    return symo

#symo = Symoro()
#robo = Robot.RX90()
#jacobian_determinant(robo, 6, range(6), range(6))
##print _jac(robo, symo, 2, 5, 5)
##print _jac_det(robo, symo, 5)
##W = kinematic_loop_constraints(robo, symo)
##print W[0]
##print W[1]
##speeds_accelerations(robo, symo)
##print _jac_inv(Robot.RX90(), symo, 2, 5, 5)
#
#def b():
#    symo = Symoro()
#    print _jac_inv(Robot.RX90(), symo, 2,5,5)
###from timeit import timeit
####print timeit(a, number=10)
####print timeit(b, number=10)
###
#import profile
##
##profile.run('b()', sort = 'cumtime')
##profile.run('b()')
#from timeit import timeit
#print timeit(b, number=1)
