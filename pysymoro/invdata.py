# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package provides symbolic
solutions for inverse geometric problem using Pieper Method.
"""


from sympy import var, sin, cos, atan2, atan, sqrt, pi
from sympy import Matrix, trigsimp, zeros
from numpy import dot, array

from pysymoro.geometry import DGM, _rot, Transform
from symoroutils import tools

EMPTY = var("EMPTY")

T_GENERAL = Matrix([var('sx,nx,ax,px'),
                    var('sy,ny,ay,py'),
                    var('sz,nz,az,pz'),
                    [0, 0, 0, 1]])


def sin_alphaj_eq_0(robo, symo, X_joint, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: Sin(alphaj) == 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.31-1.34
    # i and j revolute joints
    [i,j,k] = X_joint
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i and j are both revolute")
    symo.write_line("# Case: sin(alpha{0}) = 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        qk = robo.theta[k]
        eq_type = 3
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        el1 = cos(robo.alpha[j])*tc[z]
        el2 = cos(robo.alpha[j])*ts[z]
        el3 = cos(robo.alpha[j])*t1[z] - G[z]
        coef = [el1,el2,el3]
    else:
        qk = robo.r[k]
        eq_type = 1
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        el1 = cos(robo.alpha[j])*tr[z]
        el2 = cos(robo.alpha[j])*t2[z] - G[z]
        coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if robo.sigma[k] == 0:
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qj
    qj = robo.theta[j]
    eq_type = 3
    el1 = 2*robo.d[j]*F[x]
    el2 = -2*robo.d[j]*F[y]
    el3 = F[x]**2 + F[y]**2 + robo.d[j]**2 - (G[x]**2 + G[y]**2)
    coef = [el1,el2,el3]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    # Solve qi
    qi = robo.theta[i]
    eq_type = 4
    el1 = G[x]
    el2 = G[y]
    el3 = -robo.d[j] -cos(qj)*F[x] + sin(qj)*F[y]
    el4 = G[y]
    el5 = -G[x]
    el6 = -F[x]*sin(qj)*cos(robo.alpha[j]) - F[y]*cos(qj)*cos(robo.alpha[j])
    coef = [el1,el2,el3,el4,el5,el6]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    return

def dj_eq_0(robo, symo, X_joint, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: dj == 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.35-1.38
    # i and j revolute joints
    [i,j,k] = X_joint
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i and j are both revolute")
    symo.write_line("# Case: d{0} = 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        qk = robo.theta[k]
        eq_type = 3
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        a = 2*(tc[x]*t1[x] + tc[y]*t1[y] + tc[z]*t1[z])
        b = 2*(ts[x]*t1[x] + ts[y]*t1[y] + ts[z]*t1[z])
        c = t1[x]**2 + t1[y]**2 + t1[z]**2 + tc[x]**2 + tc[y]**2 + tc[z]**2 - (G[x]**2 + G[y]**2 + G[z]**2)
        coef = [a,b,c]
    else:
        qk = robo.r[k]
        eq_type = 2
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        a = 1
        b = 2*(tr[x]*t2[x] + tr[y]*t2[y] + tr[z]*t2[z])
        c = t2[x]**2 + t2[y]**2 + t2[z]**2 - (G[x]**2 + G[y]**2 + G[z]**2)
        coef = [a,b,c]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if robo.sigma[k] == 0:
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qj
    qj = robo.theta[j]
    eq_type = 3
    el1 = sin(robo.alpha[j])*F[y]
    el2 = sin(robo.alpha[j])*F[x]
    el3 = cos(robo.alpha[j])*F[z] - G[z]
    coef = [el1,el2,el3]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    # Solve qi
    qi = robo.theta[i]
    eq_type = 4
    el1 = G[x]
    el2 = G[y]
    el3 = -cos(qj)*F[x] + sin(qj)*F[y]
    el4 = G[y]
    el5 = -G[x]
    el6 = sin(robo.alpha[j])*F[z] - F[x]*sin(qj)*cos(robo.alpha[j]) - F[y]*cos(qj)*cos(robo.alpha[j])
    coef = [el1,el2,el3,el4,el5,el6]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    return

def dj_and_sin_alpha_dif_0(robo, symo, X_joint, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: Sin(alphaj) != 0 and dj != 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.39-1.48
    # i and j revolute joints
    [i,j,k] = X_joint
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i and j are both revolute")
    symo.write_line("# Case: d{0} != 0".format(j) + "and sin(alpha{0}) != 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        qk = robo.theta[k]
        eq_type = 6
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        K0 = t1[x]**2 + t1[y]**2 + t1[z]**2 + tc[x]**2 + tc[y]**2 + tc[z]**2 - (G[x]**2 + G[y]**2 + G[z]**2) - robo.d[j]**2
        K0 = symo.replace(trigsimp(K0), 'K0', qk)

        K1 = 2*(ts[x]*t1[x] + ts[y]*t1[y] + ts[z]*t1[z])
        K1 = symo.replace(trigsimp(K1), 'K1', qk)

        K2 = 2*(tc[x]*t1[x] + tc[y]*t1[y] + tc[z]*t1[z])
        K2 = symo.replace(trigsimp(K2), 'K2', qk)

        K3 = t1[z] - G[z]*cos(robo.alpha[j])
        K3 = symo.replace(trigsimp(K3), 'K3', qk)

        K4 = ts[z]
        K4 = symo.replace(trigsimp(K4), 'K4', qk)

        K5 = tc[z]
        K5 = symo.replace(trigsimp(K5), 'K5', qk)

        K6 = G[x]**2 + G[y]**2
        K6 = symo.replace(trigsimp(K6), 'K6', qk)

        a4 = (sin(robo.alpha[j])**2)*(K1**2 - K2**2) + 4*(robo.d[j]**2)*(K4**2 - K5**2)
        a3 = 2*(sin(robo.alpha[j])**2)*K1*K2 + 8*(robo.d[j]**2)*K4*K5
        a2 = 2*(sin(robo.alpha[j])**2)*K2*K0 + 8*(robo.d[j]**2)*K5*K3
        a1 = 2*(sin(robo.alpha[j])**2)*K1*K0 + 8*(robo.d[j]**2)*K4*K3
        a0 = (sin(robo.alpha[j])**2)*(K0**2 + K2**2) + 4*(robo.d[j]**2)*(K3**2 + K5**2) - 4*(robo.d[j]**2)*(sin(robo.alpha[j])**2)*K6
        coef = [a0,a1,a2,a3,a4]
    else:
        qk = robo.r[k]
        eq_type = 5
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        N0 = t2[x]**2 + t2[y]**2 + t2[z]**2 - (G[x]**2 + G[y]**2 + G[z]**2) - robo.d[j]**2
        N0 = symo.replace(trigsimp(N0), 'N0', qk)

        N1 = 2*(tr[x]*t2[x] + tr[y]*t2[y] + tr[z]*t2[z])
        N1 = symo.replace(trigsimp(N1), 'N1', qk)

        N2 = cos(robo.alpha[j])*G[z] + t2[z]
        N2 = symo.replace(trigsimp(N2), 'N2', qk)

        N3 = tr[z]
        N3 = symo.replace(trigsimp(N3), 'N3', qk)

        N4 = G[x]**2 + G[y]**2
        N4 = symo.replace(trigsimp(N4), 'N4', qk)

        a4 = sin(robo.alpha[j])**2
        a3 = 2*(sin(robo.alpha[j])**2)*N1
        a2 = (sin(robo.alpha[j])**2)*(2*N0 + N1**2) + 4*(robo.d[j]**2)*(N3**2)
        a1 = 2*(sin(robo.alpha[j])**2)*N1*N0 + 8*(robo.d[j]**2)*N3*N2
        a0 = (sin(robo.alpha[j])**2)*(N0**2) + 4*(robo.d[j]**2)*(N2**2) - 4*(robo.d[j]**2)*(sin(robo.alpha[j])**2)*N4
        coef = [a0,a1,a2,a3,a4]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if robo.sigma[k] == 0:
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qj
    qj = robo.theta[j]
    eq_type = 4
    el1 = sin(robo.alpha[j])*F[y]
    el2 = sin(robo.alpha[j])*F[x]
    el3 = cos(robo.alpha[j])*F[z] - G[z]
    el4 = 2*robo.d[j]*F[x]
    el5 = -2*robo.d[j]*F[y]
    el6 = robo.d[j]**2 + F[x]**2 + F[y]**2 + F[z]**2 - (G[x]**2 + G[y]**2 + G[z]**2)
    coef = [el1,el2,el3,el4,el5,el6]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    # Solve qi
    qi = robo.theta[i]
    eq_type = 4
    el1 = G[x]
    el2 = G[y]
    el3 = -cos(qj)*F[x] + sin(qj)*F[y] - robo.d[j]
    el4 = G[y]
    el5 = -G[x]
    el6 = sin(robo.alpha[j])*F[z] - F[x]*sin(qj)*cos(robo.alpha[j]) - F[y]*cos(qj)*cos(robo.alpha[j])
    coef = [el1,el2,el3,el4,el5,el6]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    return

def cos_alpha_equal_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: Cos(alphaj) == 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.50-1.51
    # i is revolute and j is prismatic
    [i,j,k] = X_joints
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i is revolute and j is prismatic")
    symo.write_line("# Case: cos(alpha{0}) = 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        # Type 3 in rk
        qk = robo.theta[k]
        eq_type = 3
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        el1 = tc[y]
        el2 = ts[y]
        el3 = t1[y] - sin(robo.alpha[j])*G[z]
        coef = [el1,el2,el3]
    else:
        # Type 1 in rk
        qk = robo.r[k]
        eq_type = 1
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        el1 = tr[y]
        el2 = t2[y] - sin(robo.alpha[j])*G[z]
        coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if robo.sigma[k] == 0:
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qi
    qi = robo.theta[i]
    eq_type = 3
    el1 = G[x]
    el2 = G[y]
    el3 = -(F[x] + robo.d[j])
    coef = [el1,el2,el3]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    # Solve qj
    qj = robo.r[j]
    eq_type = 1
    el1 = 1
    el2 = F[z] + sin(robo.alpha[j])*(-sin(qi)*G[x] + cos(qi)*G[y])
    coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    return

def cos_alpha_dif_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: Cos(alphaj) != 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.52-1.54
    # i is revolute and j is prismatic
    [i,j,k] = X_joints
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i is revolute and j is prismatic")
    symo.write_line("# Case: cos(alpha{0}) != 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        qk = robo.theta[k]
        # eq 1.53
        eq_type = 6
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        el1 = (cos(robo.alpha[j])**2)*((t1[x] + robo.d[j])**2 + tc[x]**2 - (G[x]**2 + G[y]**2)) + (t1[y] - sin(robo.alpha[j])*G[z])**2 + tc[y]**2
        el2 = 2*(cos(robo.alpha[j])**2)*(t1[x] + robo.d[j])*ts[x] + 2*(t1[y] - sin(robo.alpha[j])*G[z])*ts[y]
        el3 = 2*(cos(robo.alpha[j])**2)*(t1[x] + robo.d[j])*tc[x] + 2*(t1[y] - sin(robo.alpha[j])*G[z])*tc[y]
        el4 = 2*(cos(robo.alpha[j])**2)*tc[x]*ts[x] + 2*ts[x]*tc[y]
        el5 = (cos(robo.alpha[j])**2)*(ts[x]**2 - tc[x]**2) + ts[y]**2 - tc[y]**2
        coef = [el1,el2,el3,el4,el5]
    else:
        qk = robo.r[k]
        # eq 1.54
        eq_type = 2
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        el1 = (cos(robo.alpha[j])**2)*tr[x]**2 + tr[y]**2
        el2 = 2*(cos(robo.alpha[j])**2)*(t2[x] + robo.d[j])*tr[x] + 2*(t2[y] - sin(robo.alpha[j])*G[z])*tr[y]
        el3 = (cos(robo.alpha[j])**2)*((t2[x] + robo.d[j])**2 - (G[x]**2 + G[y]**2)) + (t2[y] - sin(robo.alpha[j])*G[z])**2
        coef = [el1,el2,el3]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if robo.sigma[k] == 0:
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qi
    qi = robo.theta[i]
    eq_type = 4
    el1 = G[x]
    el2 = G[y]
    el3 = -(F[x] + robo.d[j])
    el4 = cos(robo.alpha[j])*G[y]
    el5 = -cos(robo.alpha[j])*G[x]
    el6 = -F[y] + sin(robo.alpha[j])*G[z]
    coef = [el1,el2,el3,el4,el5,el6]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    # Solve qj
    qj = robo.r[j]
    eq_type = 1
    el1 = 1
    el2 = F[z] + sin(robo.alpha[j])*(-sin(qi)*G[x] + cos(qi)*G[y]) - cos(robo.alpha[j])*G[z]
    coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    return

def cos_alpha_equal_zero(robo, symo, X_joints, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: Cos(alphaj) == 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.56-1.57
    # i is prismatic and j is revolute
    [i,j,k] = X_joints
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i is prismatic and j is revolute")
    symo.write_line("# Case: cos(alpha{0}) = 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        qk = robo.theta[k]
        eq_type = 3
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        el1 = sin(robo.alpha[j])*tc[z]
        el2 = sin(robo.alpha[j])*ts[z]
        el3 = sin(robo.alpha[j])*t1[z] + G[y]
        coef = [el1,el2,el3]
    else:
        qk = robo.r[k]
        eq_type = 1
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        el1 = sin(robo.alpha[j])*tr[z]
        el2 = sin(robo.alpha[j])*t2[z] + G[y]
        coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if (robo.sigma[k] == 0):
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qj
    qj = robo.theta[j]
    eq_type = 3
    el1 = F[x]
    el2 = -F[y]
    el3 = -G[x] + robo.d[j]
    coef = [el1,el2,el3]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    # Solve qi
    qi = robo.r[i]
    eq_type = 1
    el1 = 1
    el2 = -G[z] + sin(robo.alpha[j])*(sin(qj)*F[x] + cos(qj)*F[y])
    coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    return

def cos_alpha_dif_zero(robo, symo, X_joints, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.
    (Case: Cos(alphaj) != 0)
    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.58-1.61
    # i is prismatic and j is revolute
    [i,j,k] = X_joints
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i is prismatic and j is revolute")
    symo.write_line("# Case: cos(alpha{0}) != 0".format(j) + "\r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        eq_type = 6
        qk = robo.theta[k]
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        K3 = sin(robo.alpha[j])*t1[z] + G[y]
        K3 = symo.replace(trigsimp(K3), 'K3', qk)

        K4 = sin(robo.alpha[j])*ts[z]
        K4 = symo.replace(trigsimp(K4), 'K4', qk)

        K5 = sin(robo.alpha[j])*tc[z]
        K5 = symo.replace(trigsimp(K5), 'K5', qk)

        K6 = t1[x]**2 + t1[y]**2
        K6 = symo.replace(trigsimp(K6), 'K6', qk)

        K7 = 2*(ts[x]*t1[x] + ts[y]*t1[y])
        K7 = symo.replace(trigsimp(K7), 'K7', qk)

        K8 = 2*(tc[x]*t1[x] + tc[y]*t1[y])
        K8 = symo.replace(trigsimp(K8), 'K8', qk)

        K9 = 2*(tc[x]*ts[x] + tc[y]*ts[y])
        K9 = symo.replace(trigsimp(K9), 'K9', qk)

        K10 = ts[x]**2 + ts[y]**2
        K10 = symo.replace(trigsimp(K10), 'K10', qk)

        K11 = tc[x]**2 + tc[y]**2
        K11 = symo.replace(trigsimp(K11), 'K11', qk)

        el1 = K3**2 + K5**2 + (cos(robo.alpha[j])**2)*((G[x] - robo.d[j])**2 - (K6 + K11))
        el2 = 2*K4*K3 - (cos(robo.alpha[j])**2)*K7
        el3 = 2*K5*K3 - (cos(robo.alpha[j])**2)*K8
        el4 = 2*K4*K5 - (cos(robo.alpha[j])**2)*K9
        el5 = K4**2 - K5**2 - (cos(robo.alpha[j])**2)*(K10 - K11)
        coef = [el1,el2,el3,el4,el5]
    else:
        qk = robo.r[k]
        eq_type = 2
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        el1 = tr[z]**2 - cos(robo.alpha[j])**2
        el2 = -2*(cos(robo.alpha[j])**2)*(tr[x]*t2[x] + tr[y]*t2[y]) + 2*sin(robo.alpha[j])*tr[z]*(sin(robo.alpha[j])*t2[z] + G[y])
        el3 = (cos(robo.alpha[j])**2)*((G[x] - robo.d[j])**2 - t2[x]**2 - t2[y]**2) + (sin(robo.alpha[j])*t2[z] + G[y])**2
        coef = [el1,el2,el3]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if robo.sigma[k] == 0:
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qj
    qj = robo.theta[j]
    eq_type = 4
    el1 = F[x]
    el2 = -F[y]
    el3 = -G[x] + robo.d[j]
    el4 = cos(robo.alpha[j])*F[y]
    el5 = cos(robo.alpha[j])*F[x]
    el6 = -(G[y] + sin(robo.alpha[j])*F[z])
    coef = [el1,el2,el3,el4,el5,el6]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    # Solve qi
    qi = robo.r[i]
    eq_type = 1
    el1 = -1
    el2 = G[z] -cos(robo.alpha[j])*F[z] - sin(robo.alpha[j])*(sin(qj)*F[x] + cos(qj)*F[y])
    coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    return

def ij_prismatic(robo, symo, X_joint, tc, ts, tr, t0, G, offset):
    """
    Function that finds the symbolic solutions of the X joints.

    Parameters:
    ===========
    1) X_joints:  Type of the X joints
    2) tc, ts, tr, t0: Coefficients of eq. 1.25
    3) G: G = [Gx; Gy; Gz] -> Vector with constant values
    """
    # eq 1.62-1.63 (Corrected version)
    # i and j are prismatic joints
    [i,j,k] = X_joint
    [x,y,z] = [0,1,2]
    [offseti, offsetj, offsetk] = offset
    symo.write_line("# X joints i and j are both prismatic \r\n")

    # Solve qk
    if robo.sigma[k] == 0:
        qk = robo.theta[k]
        eq_type = 3
        t1 = tr*robo.r[k] + t0
        t1 = symo.replace(trigsimp(t1), 't1', qk)

        el1 = tc[x]
        el2 = ts[x]
        el3 = t1[x] - G[x] + robo.d[j]
        coef = [el1,el2,el3]
    else:
        qk = robo.r[k]
        eq_type = 1
        t2 = tc*cos(robo.theta[k]) + ts*sin(robo.theta[k]) + t0
        t2 = symo.replace(trigsimp(t2), 't2', qk)

        el1 = tr[x]
        el2 = t2[x] - G[x] + robo.d[j]
        coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qk,offsetk)

    if (robo.sigma[k] == 0):
        F = tc*cos(qk) + ts*sin(qk) + t1
    else:
        F = tr*qk + t2
    F = symo.replace(trigsimp(F), 'F', qk)

    # Solve qj
    qj = robo.r[j]
    eq_type = 1
    el1 = -sin(robo.alpha[j])
    el2 = cos(robo.alpha[j])*F[y] - sin(robo.alpha[j])*F[z] - G[y]
    coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qj,offsetj)

    # Solve qi
    qi = robo.r[i]
    eq_type = 1
    el1 = 1
    el2 = sin(robo.alpha[j])*F[y] + cos(robo.alpha[j])*F[z] - G[z] + cos(robo.alpha[j])*qj
    coef = [el1,el2]
    _equation_solve(symo,coef,eq_type,qi,offseti)

    return

def _equation_solve(symo, coef, eq_type, unknown, offset):
    """
    Function that solves the possible type of equations that we need to solve.
    (One extra type is the system of equations for the prismatic case -> type 0)

    Parameters:
    ===========
    1) Symo instance
    2) coef: Vector containing the coefficients of each equation
    3) eq_type: Type of equation for the unknown (System of equation is assigned as 0 type)
    4) unknown: The unknown(s) of this equation system
    """

    symo.write_line("\r\n\r\n")
    if eq_type == 0:
        """
        System of equations for the three prismatic joints case
        """
        [ri,rj,rk] = unknown
        [a1,a2,ct] = [0,1,2]
        [x,y,z] = [0,1,2]

        symo.write_line("# System of equations in " + str(ri) + ", " + str(rj) + " and " + str(rk) + ":")
        symo.write_line("#=================================================\r\n")

        a1x = symo.replace(coef[a1,x], 'a1x', rk)
        a1y = symo.replace(coef[a1,y], 'a1y', rk)
        a1z = symo.replace(coef[a1,z], 'a1z', rk)
        a2x = symo.replace(coef[a2,x], 'a2x', ri)
        a2y = symo.replace(coef[a2,y], 'a2y', ri)
        a2z = symo.replace(coef[a2,z], 'a2z', ri)
        ct1 = symo.replace(coef[ct,x], 'ct1')
        ct2 = symo.replace(coef[ct,y], 'ct2')
        ct3 = symo.replace(coef[ct,z], 'ct3')

        expr1 = array(a1x)*rk + array(a2x)*ri
        expr2 = array(a1y)*rk + array(a2y)*ri
        expr3 = array(a1z)*rk + array(a2z)*ri + rj
        symo.write_line("\r\n# Equations: ")
        symo.write_line("#========== \r\n")
        symo.write_line("# {0}".format(expr1) + " = " + "{0}".format(-array(ct1)))
        symo.write_line("# {0}".format(expr2) + " = " + "{0}".format(-array(ct2)))
        symo.write_line("# {0}".format(expr3) + " = " + "{0}".format(-array(ct3)) + "\r\n")
        symo.write_line("# Solution: ")
        symo.write_line("#=========\r\n")
        A = Matrix( [[a1x,a2x,0],[a1y,a2y,0],[a1z,a2z,1]] )
        B = -Matrix( [ct1,ct2,ct3] )
        C = Matrix( [offset[0], offset[1], offset[2]] )
        r = A.inv()*B - C
        symo.write_line("\r\n\r\n# The solutions of the lengths for the given linear system, are:")
        symo.write_line("#==============================================================\r\n\r\n")
        symo.write_line("# {0}".format(rk) + " = " + "{0}".format(r[0]))
        symo.write_line("# {0}".format(ri) + " = " + "{0}".format(r[1]))
        symo.write_line("# {0}".format(rj) + " = " + "{0}".format(r[2]))
        symo.write_line("\r\n")

    elif eq_type == 1:
        """Solution for the system:
        a*r + b = 0
        """
        r = unknown
        symo.write_line("# Equation in {0}: ".format(r))
        symo.write_line("#=================")
        symo.write_line("# Type {0}".format(eq_type))
        symo.write_line("# a*{0} + b = 0".format(r) + "\r\n")
        a = symo.replace(coef[0], "a", r)
        b = symo.replace(-coef[1], "b", r)
        expr = a*r
        symo.write_line("\r\n# Equation: ")
        symo.write_line("#========== \r\n")
        symo.write_line("# {0}".format(expr) + " = " + "{0}".format(b) + "\r\n")
        symo.write_line("# Solution: ")
        symo.write_line("#==============\r\n")
        symo.add_to_dict(r, b/a - offset)
        symo.write_line("\r\n")

    elif eq_type == 2:
        """Solution for the system:
        a*r^2 + b*r + c = 0
        """
        r = unknown
        symo.write_line("# Equation in {0}: ".format(r))
        symo.write_line("#=================")
        EPS = var('EPS'+str(r))
        symo.write_line("# Type {0}".format(eq_type))
        symo.write_line("# a*{0}**2 + b*{0} + c = 0".format(r) + "\r\n")
        a = symo.replace(coef[0],"a",r)
        b = symo.replace(coef[1],"b",r)
        c = symo.replace(coef[2],"c",r)
        expr = a*r**2 + b*r + c
        symo.write_line("\r\n# Equation: ")
        symo.write_line("#========== \r\n")
        symo.write_line("# {0} = 0".format(expr) + "\r\n")
        symo.add_to_dict(EPS, (tools.ONE, - tools.ONE))
        Delta = symo.replace(b**2 - 4*a*c, 'Delta', r)
        symo.write_line("# Solution: ")
        symo.write_line("#==============\r\n")
        symo.add_to_dict(r, (-b + EPS*sqrt(Delta))/2*a - offset)
        symo.write_line("\r\n")

    elif eq_type == 3:
        """Solution for the system:
        a*C(th) + b*S(th) + c = 0
        """
        th = unknown
        symo.write_line("# Equation in {0}: ".format(th))
        symo.write_line("#=================")
        EPS = var('EPS'+str(th))
        symo.write_line("# Type {0}".format(eq_type))
        symo.write_line("# a*C({0}) + b*S({0}) + c = 0".format(th) + "\r\n")
        a = symo.replace(coef[0], "a", th)
        b = symo.replace(coef[1], "b", th)
        c = symo.replace(-coef[2], "c", th)
        expr = a*cos(th) + b*sin(th)
        symo.write_line("\r\n# Equation: ")
        symo.write_line("#========== \r\n")
        symo.write_line("# {0}".format(expr) + " = " + " {0}".format(c) + "\r\n")
        symo.write_line("# Solution: ")
        symo.write_line("#==============\r\n")

        # Type 3
        # case 1
        if a == tools.ZERO and b != tools.ZERO:
            S = symo.replace(c/b, 'S', th)
            symo.add_to_dict(EPS, (tools.ONE, - tools.ONE))
            symo.add_to_dict(th, atan2(S, EPS*sqrt(1-S**2)) - offset )
        # case 2
        elif a != tools.ZERO and b == tools.ZERO:
            C = symo.replace(c/a, 'C', th)
            symo.add_to_dict(EPS, (tools.ONE, - tools.ONE))
            symo.add_to_dict(th, atan2(EPS*sqrt(1-C**2), C) - offset )
        # case 3
        elif c == tools.ZERO:
            symo.add_to_dict(EPS, (tools.ONE, tools.ZERO))
            symo.add_to_dict(th, atan2(-a, b) + EPS*pi - offset )
        # case 4
        else:
            B = symo.replace(a**2 + b**2, 'B', th)
            D = symo.replace(B - c**2, 'D', th)
            symo.add_to_dict(EPS, (tools.ONE, - tools.ONE))
            S = symo.replace((b*c + EPS*a*sqrt(D))/B, 'S', th)
            C = symo.replace((a*c - EPS*b*sqrt(D))/B, 'C', th)
            symo.add_to_dict(th, atan2(S, C) - offset )
        symo.write_line("\r\n")

    elif eq_type == 4:
        """Solution for the system:
        a*C(th) + b*S(th) + c = 0
        a'*C(th) + b'*S(th) + c' = 0
        """
        th = unknown
        symo.write_line("# Equation in {0}: ".format(th))
        symo.write_line("#=================")
        symo.write_line("# Type {0}".format(eq_type))
        symo.write_line("# a*C({0}) + b*S({0}) + c = 0".format(th))
        symo.write_line("# a'*C({0}) + b'*S({0}) + c' = 0".format(th) + "\r\n")
        a = symo.replace(coef[0], 'a', th)
        b = symo.replace(coef[1], 'b', th)
        c = symo.replace(coef[2], 'c', th)
        ap = symo.replace(coef[3], 'ap', th)
        bp = symo.replace(coef[4], 'bp', th)
        cp = symo.replace(coef[5], 'cp', th)
        expr1 = a*cos(th) + b*sin(th)
        expr2 = ap*cos(th) + bp*sin(th)
        symo.write_line("\r\n# Equations:".format(th))
        symo.write_line("#===============")
        symo.write_line("# {0}".format(expr1) + " = " + "{0}".format(-c) )
        symo.write_line("# {0}".format(expr2) + " = " + "{0}".format(-cp) + "\r\n")
        symo.write_line("# Solution: ")
        symo.write_line("#==============\r\n")
        if b == tools.ZERO and ap == tools.ZERO:
            symo.add_to_dict(th, atan2(-cp/bp, -c/a) - offset )
        elif bp == tools.ZERO and a == tools.ZERO:
            symo.add_to_dict(th, atan2(-c/b, -cp/ap) - offset )
        else:
            D = symo.replace(b*ap - bp*a, 'D', th)
            C = symo.replace(-(cp*b - c*bp)/D, 'C', th)
            S = symo.replace(-(c*ap - cp*a)/D, 'S', th)
            symo.add_to_dict(th, atan2(S, C) - offset)
        symo.write_line("\r\n")

    elif eq_type == 5:
        """Solution for the system:
        a4*r^4 + a3*r^3 + a2*r^2 + a1*r + a0 = 0
        """
        r = unknown
        EPS = var('EPS'+str(r))
        symo.write_line("\r\n# Equation in {0}: ".format(r))
        symo.write_line("#=================")
        symo.write_line("# Type {0}".format(eq_type))
        symo.write_line("# a4*{0}**4 + a3*{0}**3 + a2*{0}**2 + a1*{0} + a0 = 0".format(r) + "\r\n")
        a4 = symo.replace(coef[4], "a4", r)
        a3 = symo.replace(coef[3], "a3", r)
        a2 = symo.replace(coef[2], "a2", r)
        a1 = symo.replace(coef[1], "a1", r)
        a0 = symo.replace(coef[0], "a0", r)
        expr = a4*r**4 + a3*r**3 + a2*r**2 + a1*r + a0
        symo.write_line("\r\n# Equations:".format(r))
        symo.write_line("#===============")
        symo.write_line("# {0} = 0".format(expr) )
        symo.write_line("\r\n# Solution: ")
        symo.write_line("#==============\r\n")
        delta0 = symo.replace(trigsimp( a2**2 - 3*a3*a1 + 12*a4*a0 ), "delta0", r)
        delta1 = symo.replace(trigsimp( 2*a2**3 - 9*a3*a2*a1 + 27*(a3**2)*a0 + 27*a4*a1**2 - 72*a4*a2*a0 ), "delta1", r)
        Q = symo.replace(trigsimp( ((delta1 + sqrt(delta1**2 - 4*delta0**3))/2)**(1/3) ), "Q", r)
        p = symo.replace(trigsimp( (8*a4*a2 - 3*a3**2)/(8*a4**2) ), "p", r)
        q = symo.replace(trigsimp( (a3**3 - 4*a4*a3*a2 + 8*(a4**2)*a1)/(8*a4**3) ), "q", r)
        S = symo.replace(trigsimp( sqrt(-(2*p)/3 + ((Q**2 + delta0)/(3*a4*Q)))/2 ), "S", r)
        symo.add_to_dict(EPS, (tools.ONE, - tools.ONE))
        sol12 = -a3/(4*a4) - S + EPS*sqrt(-4*S**2 - 2*p + q/S)/2 - offset
        sol34 = -a3/(4*a4) + S + EPS*sqrt(-4*S**2 - 2*p - q/S)/2 - offset
        symo.write_line("\r\n# Solutions 1 and 2 are: ")
        r_12 = var(str(r)+'_12')
        symo.add_to_dict(r_12, sol12)
        symo.write_line("\r\n\r\n# Solutions 3 and 4 are: \r\n")
        r_34 = var(str(r)+'_34')
        symo.add_to_dict(r_34, sol34)

    elif eq_type == 6:
        """Solution for the system:
        a4*S(th)^2 + a3*C(th)*S(th) + a2*C(th) + a1*S(th) + a0 = 0
        Equation in type 5:
        (a0 - a2)*t^4 + 2*(a1 - a3)*t^3 + 2*(a0 + 2*a4)*t^2 + 2*(a1 + a3)*t + (a0 + a2) = 0 , with t = tan(th/2)
        """
        th = unknown
        t = var('t')
        EPS = var('EPS'+str(th))
        symo.write_line("# Equation in {0}: ".format(th))
        symo.write_line("#=================")
        symo.write_line("# Type {0}".format(eq_type))
        symo.write_line("# a4*S({0})**2 + a3*CS({0}) + a2*C({0}) + a1*S({0}) + a0 = 0".format(th) + "\r\n")
        a4 = symo.replace(coef[4], "a4", th)
        a3 = symo.replace(coef[3], "a3", th)
        a2 = symo.replace(coef[2], "a2", th)
        a1 = symo.replace(coef[1], "a1", th)
        a0 = symo.replace(coef[0], "a0", th)
        expr1 = a4*sin(th)**2 + a3*cos(th)*sin(th) + a2*cos(th) + a1*sin(th) + a0
        A4 = symo.replace(trigsimp(a0 - a2), "A4", th)
        A3 = symo.replace(trigsimp(2*(a1 - a3)), "A3", th)
        A2 = symo.replace(trigsimp(2*(a0 + 2*a4)), "A2", th)
        A1 = symo.replace(trigsimp(2*(a1 + a3)), "A1", th)
        A0 = symo.replace(trigsimp(a0 + a2), "A0", th)
        expr2 = A4*t**4 + A3*t**3 + A2*t**2 + A1*t + A0
        symo.write_line("\r\n# Equation:".format(th))
        symo.write_line("#===============")
        symo.write_line("# {0} = 0".format(expr1) )
        symo.write_line("\r\n# New Equation:")
        symo.write_line("#=================")
        symo.write_line("# {0} = 0".format(expr2) + "  (with {0}".format(t) + " = " + "tan({0}/2) )".format(th) + "\r\n\r\n")
        symo.write_line("\r\n# Solution: ")
        symo.write_line("#==============\r\n")
        delta0 = symo.replace(trigsimp( A2**2 - 3*A3*A1 + 12*A4*A0 ), "delta0", th)
        delta1 = symo.replace(trigsimp( 2*A2**3 - 9*A3*A2*A1 + 27*(A3**2)*A0 + 27*A4*A1**2 - 72*A4*A2*A0 ), "delta1", th)
        Q = symo.replace(trigsimp( ((delta1 + sqrt(delta1**2 - 4*delta0**3))/2)**(1/3) ), "Q", th)
        p = symo.replace(trigsimp( (8*A4*A2 - 3*A3**2)/(8*A4**2) ), "p", th)
        q = symo.replace(trigsimp( (A3**3 - 4*A4*A3*A2 + 8*(A4**2)*A1)/(8*A4**3) ), "q", th)
        S = symo.replace(trigsimp( sqrt(-(2*p)/3 + ((Q**2 + delta0)/(3*A4*Q)))/2 ), "S", th)
        symo.add_to_dict(EPS, (tools.ONE, - tools.ONE))
        sol12 = -A3/(4*A4) - S + EPS*sqrt(-4*S**2 - 2*p + q/S)/2
        sol34 = -A3/(4*A4) + S + EPS*sqrt(-4*S**2 - 2*p - q/S)/2
        symo.write_line("\r\n# Solutions 1 and 2 are: ")
        t_12 = symo.replace(sol12, "t_12", th)
        symo.add_to_dict("{0}_12".format(th), 2*atan(t_12) - offset )
        symo.write_line("\r\n# Solutions 3 and 4 are: \r\n")
        t_34 = symo.replace(sol34, "t_34", th)
        symo.add_to_dict("{0}_34".format(th), 2*atan(t_34) - offset )

    symo.write_line("--------------------------------------------------------------------------------------------")
    return

def solve_position(robo, symo, com_key, X_joints, fc, fs, fr, f0, g):
    """
    Function that solves the position equation for the four spherical cases.

    Parameters:
    ===========
    1) com_key: X joints combination
    2) X_joints: Type of the X joints
    3) fc, fs, fr, f0: Coefficients of equation (1.23.b)
    4) g: g = [gx; gy; gz] -> Vector containing constant values
    """
    [i,j,k] = X_joints
    [x,y,z] = [0,1,2]           # Book convention indexes
    offset = zeros(1,len(X_joints))

    if robo.sigma[k] == 0:
        robo.r[k] = robo.r[k] + robo.b[robo.ant[k]]
        offset[2] = robo.gamma[robo.ant[k]]
    elif robo.sigma[k] == 1:
        robo.theta[k] = robo.theta[k] + robo.gamma[robo.ant[k]]
        offset[2] = robo.b[robo.ant[k]]

    if robo.sigma[j] == 0:   # If the joint j is revolute
        robo.r[j] = robo.r[j] + robo.b[robo.ant[j]]
        offset[1] = robo.gamma[robo.ant[j]]
        T = Transform.create(axis=z, th=0, p=robo.r[j])
    elif robo.sigma[j] == 1:                         # If the joint j is prismatic
        robo.theta[j] = robo.theta[j] + robo.gamma[robo.ant[j]]
        offset[1] = robo.b[robo.ant[j]]
        T = Transform.create(axis=z, th=robo.theta[j], p=0)

    tc = T*fc
    tc = symo.replace(trigsimp(tc), 'tc', k)
    ts = T*fs
    ts = symo.replace(trigsimp(ts), 'ts', k)
    tr = T*fr
    tr = symo.replace(trigsimp(tr), 'tr', k)
    t0 = T*f0
    t0 = symo.replace(trigsimp(t0), 't0', k)

    if robo.sigma[i] == 0:   # If joint i is revolute
        robo.r[i] = robo.r[i] + robo.b[robo.ant[i]]
        offset[0] = robo.gamma[robo.ant[i]]
        T = Transform.create(axis=z, th=0, p=-robo.r[i])
    else:                         # If the joint i is prismatic
        robo.theta[i] = robo.theta[i] + robo.gamma[robo.ant[i]]
        offset[0] = robo.b[robo.ant[i]]
        T = Transform.create(axis=z, th=-robo.theta[i], p=0)

    G = T*g
    G = symo.replace(trigsimp(G), 'G')

    # Conditions to get the solution for its possible X joints combination
    if (com_key == 0) or (com_key == 1):                                  # X joints: RRX with X:either R or P joint
        if sin(robo.alpha[j]) == 0:                              # sin(alphaj) equal 0
            sin_alphaj_eq_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
        elif robo.d[j] == 0:                                     # dj equal 0
            dj_eq_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
        elif (sin(robo.alpha[j]) != 0) and (robo.d[j] != 0):     # sin(alphaj) and dj not equal to 0
            dj_and_sin_alpha_dif_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
    elif (com_key == 2) or (com_key == 3):                                # X joints: RPX with X:either R or P joint
        if cos(robo.alpha[j]) == 0:                              # cos(alphaj) equal 0
            cos_alpha_equal_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
        elif cos(robo.alpha[j]) != 0:                            # cos(alphaj) not equal to 0
            cos_alpha_dif_0(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
    elif (com_key == 4) or (com_key == 5):                                # X joints: PRX with X:either R or P joint
        if cos(robo.alpha[j]) == 0:                              # cos(alphaj) equal 0
            cos_alpha_equal_zero(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
        elif cos(robo.alpha[j]) != 0:                            # cos(alphaj) not equal to 0
            cos_alpha_dif_zero(robo, symo, X_joints, tc, ts, tr, t0, G, offset)
    else:                                                                 # X joints: PPX with X:either R or P joint
        ij_prismatic(robo, symo, X_joints, tc, ts, tr, t0, G, offset)

    return

def solve_orientation(robo, symo, pieper_joints):
    """
    Function that solves the orientation equation for the four spherical cases.

    Parameters:
    ===========
    1) pieper_joints: Joints that form the spherical wrist
    """
    m = pieper_joints[1]        # Center of the spherical joint
    [S,N,A] = [0,1,2]           # Book convention indexes
    [x,y,z] = [0,1,2]           # Book convention indexes

    t1 = dgm(robo, symo, m-2, 0, fast_form=True, trig_subs=True)
    t2 = dgm(robo, symo, 6, m+1, fast_form=True, trig_subs=True)

    Am2A0 = Matrix([ t1[:3,:3] ])
    A6Am1 = Matrix([ t2[:3,:3] ])

    A0 = T_GENERAL[:3,:3]

    SNA = _rot(axis=x, th=-robo.alpha[m-1])*Am2A0*A0*A6Am1
    SNA = symo.replace(trigsimp(SNA), 'SNA')

    # calculate theta(m-1) (i)
    # eq 1.68
    eq_type = 3
    offset = robo.gamma[robo.ant[m-1]]
    robo.r[m-1] = robo.r[m-1] + robo.b[robo.ant[m-1]]
    coef = [-SNA[y,A]*sin(robo.alpha[m]) , SNA[x,A]*sin(robo.alpha[m]) , SNA[z,A]*cos(robo.alpha[m])-cos(robo.alpha[m+1])]
    _equation_solve(symo, coef, eq_type, robo.theta[m-1], offset)

    # calculate theta(m) (j)
    # eq 1.72
    S1N1A1 = _rot(axis=x, th=-robo.alpha[m])*_rot(axis=z, th=-robo.theta[m-1])*SNA
    eq_type = 4
    offset = robo.gamma[robo.ant[m]]
    robo.r[m] = robo.r[m] + robo.b[robo.ant[m]]
    symo.write_line("\r\n\r\n")
    B1 = symo.replace(trigsimp(-S1N1A1[x,A]), 'B1', robo.theta[m])
    B2 = symo.replace(trigsimp(S1N1A1[y,A]), 'B2', robo.theta[m])
    coef = [0, sin(robo.alpha[m+1]), B1, sin(robo.alpha[m+1]), 0, B2]
    _equation_solve(symo, coef, eq_type, robo.theta[m], offset)

    # calculate theta(m+1) (k)
    # eq 1.73
    eq_type = 4
    offset = robo.gamma[robo.ant[m+1]]
    robo.r[m+1] = robo.r[m+1] + robo.b[robo.ant[m+1]]
    symo.write_line("\r\n\r\n")
    B1 = symo.replace(trigsimp(-S1N1A1[z,S]), 'B1', robo.theta[m+1])
    B2 = symo.replace(trigsimp(-S1N1A1[z,N]), 'B2', robo.theta[m+1])
    coef = [0, sin(robo.alpha[m+1]), B1, sin(robo.alpha[m+1]), 0, B2]
    _equation_solve(symo, coef, eq_type, robo.theta[m+1], offset)

    return


def solve_orientation_prismatic(robo, symo, X_joints):
    """
    Function that solves the orientation equation of the three prismatic joints case. (to find the three angles)

    Parameters:
    ===========
    1) X_joint: The three revolute joints for the prismatic case
    """
    [i,j,k] = X_joints                   # X joints vector
    [S,S1,S2,S3,x] = [0,0,0,0,0]        # Book convention indexes
    [N,N1,N2,N3,y] = [1,1,1,1,1]        # Book convention indexes
    [A,A1,A2,A3,z] = [2,2,2,2,2]        # Book convention indexes
    robo.theta[i] = robo.theta[i] + robo.gamma[robo.ant[i]]
    robo.theta[j] = robo.theta[j] + robo.gamma[robo.ant[j]]
    robo.theta[k] = robo.theta[k] + robo.gamma[robo.ant[k]]

    T6k = dgm(robo, symo, 6, k, fast_form=True,trig_subs=True)
    Ti0 = dgm(robo, symo, robo.ant[i], 0, fast_form=True, trig_subs=True)
    Tji = dgm(robo, symo, j-1, i, fast_form=True, trig_subs=True)
    Tjk = dgm(robo, symo, j, k-1, fast_form=True, trig_subs=True)

    S3N3A3 = Transform.create(axis=x, th=-robo.alpha[i], p=0)*Ti0*T_GENERAL*T6k   # S3N3A3 = rot(x,-alphai)*rot(z,-gami)*aiTo*SNA
    S2N2A2 = Transform.create(axis=x, th=-robo.alpha[j], p=0)*Tji                 # S2N2A2 = iTa(j)*rot(x,-alphaj)
    S1N1A1 = Tjk*Transform.create(axis=x, th=-robo.alpha[k], p=0)                 # S1N1A1 = jTa(k)*rot(x,alphak)

    S3N3A3 = Matrix([ S3N3A3[:3, :3] ])
    S2N2A2 = Matrix([ S2N2A2[:3, :3] ])
    S1N1A1 = Matrix([ S1N1A1[:3, :3] ])
    SNA = Matrix([ T_GENERAL[:3, :3] ])
    SNA = symo.replace(trigsimp(SNA), 'SNA')

    # solve thetai
    # eq 1.100 (page 49)
    eq_type = 3
    offset = robo.gamma[robo.ant[i]]
    robo.r[i] = robo.r[i] + robo.b[robo.ant[i]]
    el1 = S2N2A2[z,S2]*S3N3A3[x,A3] + S2N2A2[z,N2]*S3N3A3[y,A3]
    el2 = S2N2A2[z,S2]*S3N3A3[y,A3] - S2N2A2[z,N2]*S3N3A3[x,A3]
    el3 = S2N2A2[z,A2]*S3N3A3[z,A3] - S1N1A1[z,A1]
    coef = [el1,el2,el3]
    _equation_solve(symo, coef, eq_type, robo.theta[i], offset)

    # solve thetaj
    # eq 1.102
    eq_type = 4
    offset = robo.gamma[robo.ant[j]]
    robo.r[j] = robo.r[j] + robo.b[robo.ant[j]]
    coef = [S1N1A1[x,A1] , -S1N1A1[y,A1] , -SNA[x,A] , S1N1A1[y,A1] , S1N1A1[x,A1] , -SNA[y,A] ]
    _equation_solve(symo, coef, eq_type, robo.theta[j], offset)

    # solve thetak
    # eq 1.103
    eq_type = 4
    offset = robo.gamma[robo.ant[k]]
    robo.r[k] = robo.r[k] + robo.b[robo.ant[k]]
    coef = [S1N1A1[z,S1] , S1N1A1[z,N1] , -SNA[z,S] , S1N1A1[z,N1] , -S1N1A1[z,S1] , -SNA[z,N] ]
    _equation_solve(symo, coef, eq_type, robo.theta[k], offset)

    return


def solve_position_prismatic(robo, symo, pieper_joints):
    """
    Function that solves the position equation of the three prismatic joints case. (to find the three angles)

    Parameters:
    ===========
    1) Pieper_joints: The three prismatic joints for the prismatic case
    """
    eq_type = 0                       # Type 0: Linear system
    [i,j,k] = pieper_joints           # Pieper joints vector
    [x,y,z] = [0,1,2]                 # Book convention indexes
    robo.theta[i] = robo.theta[i] + robo.gamma[robo.ant[i]]
    robo.theta[j] = robo.theta[j] + robo.gamma[robo.ant[j]]
    robo.theta[k] = robo.theta[k] + robo.gamma[robo.ant[k]]

    T6k = dgm(robo, symo, 6, k, fast_form=True, trig_subs=True)
    Ti0 = dgm(robo, symo, robo.ant[i], 0, fast_form=True, trig_subs=True)
    Tji = dgm(robo, symo, j-1, i, fast_form=True, trig_subs=True)
    Tjk = dgm(robo, symo, j, k-1, fast_form=True, trig_subs=True)

    S3N3A3P3 = Transform.create(axis=z, th=-robo.theta[i], p=0)*Transform.create(axis=x, th=-robo.alpha[i], p=-robo.d[i])*Ti0*T_GENERAL*T6k # S3N3A3 = rot(z,-thetai)*Trans(x,-di)*rot(x,-alphai)*a(i)To*SNAP*6Tk
    S2N2A2P2 = Transform.create(axis=z, th=-robo.theta[j], p=0)*Transform.create(axis=x, th=-robo.alpha[j], p=-robo.d[j])*Tji               # S2N2A2 = rot(z,-thetaj)*Trans(x,-dj)*rot(x,-alphaj)*a(j)Ti
    S1N1A1P1 = Tjk*Transform.create(axis=x, th=robo.alpha[k], p=robo.d[k])*Transform.create(axis=z, th=robo.theta[k], p=0)                  # S1N1A1 = jTa(k)*rot(x,alphak)*Trans(x,dk)*rot(z,thetak)

    S2N2A2 = array( S2N2A2P2[:3, :3] )
    P3 = array( S3N3A3P3[:3, 3] )
    P2 = array( S2N2A2P2[:3, 3] )
    P1 = array( S1N1A1P1[:3, 3] )

    t2 = S2N2A2[:3, 2]
    P4 = P1 - dot(S2N2A2, P3) - P2
    t1 = S1N1A1P1[:3,2]
    coef = Matrix([[t1[0],t1[1],t1[2]],[t2[0],t2[1],t2[2]],[P4[0][0],P4[1][0],P4[2][0]]])
    rijk = [robo.r[i], robo.r[j], robo.r[k]]
    offset = [robo.b[robo.ant[k]], robo.b[robo.ant[j]], robo.b[robo.ant[i]]]
    _equation_solve(symo, coef, eq_type, rijk, offset)

    return
