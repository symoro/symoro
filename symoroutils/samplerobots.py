# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module contains some sample robots (parameters) that can be used
for multiple purposes.
"""


from sympy import pi, var, zeros
from sympy import Matrix

from pysymoro.robot import Robot
from symoroutils import tools


def cart_pole():
    """Generate Robot instance of classical CartPole dynamic system."""
    #TODO: bring it to the new notation with 0-frame
    robo = Robot('CartPole', 2, 2, 2, False)
    robo.ant = (-1, 0, 1)
    robo.sigma = (0, 1, 0)
    robo.alpha = (0, pi/2, pi/2)
    robo.d = (0, 0, 0)
    robo.theta = (0, pi/2, var('th2'))
    robo.r = (0, var('r1'), 0)
    robo.b = (0, 0, 0)
    robo.gamma = (0, 0, 0)
    robo.structure = tools.SIMPLE
    robo.num = range(0, 3)
    robo.Nex = [zeros(3, 1) for i in robo.num]
    robo.Fex = [zeros(3, 1) for i in robo.num]
    robo.FS = [0 for i in robo.num]
    robo.IA = [0 for i in robo.num]
    robo.FV = [var('FV{0}'.format(i)) for i in robo.num]
    robo.MS = [zeros(3, 1) for i in robo.num]
    robo.MS[1][0] = var('MX2')
    robo.M = [var('M{0}'.format(i)) for i in robo.num]
    robo.GAM = [var('GAM{0}'.format(i)) for i in robo.num]
    inertia_matrix_terms = ("XX{0}, XY{0}, XZ{0}, ") + \
        ("XY{0}, YY{0}, YZ{0}, ") + \
        ("XZ{0}, YZ{0}, ZZ{0}")
    robo.J = [
        Matrix(3, 3, var(inertia_matrix_terms.format(i))) \
        for i in robo.num
    ]
    robo.G = Matrix([0, 0, -var('G3')])
    robo.w0 = zeros(3, 1)
    robo.wdot0 = zeros(3, 1)
    robo.v0 = zeros(3, 1)
    robo.vdot0 = zeros(3, 1)
    robo.q = [0, var('r1'), var('th2')]
    robo.qdot = [0, var('r1d'), var('th2d')]
    robo.qddot = [0, var('r1dd'), var('th2dd')]
    return robo


def planar2r():
    """Generate Robot instance of 2R Planar robot"""
    robo = Robot('Planar2R', 2, 2, 2, False)
    robo.structure = tools.SIMPLE
    robo.sigma = [2, 0, 0]
    robo.mu = [0, 1, 1]
    robo.gamma = [0, 0, 0]
    robo.b = [0, 0, 0]
    robo.alpha = [0, 0, 0]
    robo.d = [0, 0, var('L1')]
    robo.theta = [0, var('q1'), var('q2')]
    robo.r = [0, 0, 0]
    robo.num = range(0, 3)
    robo.Nex = [zeros(3, 1) for i in robo.num]
    robo.Fex = [zeros(3, 1) for i in robo.num]
    robo.FS = [0 for i in robo.num]
    robo.IA = [0 for i in robo.num]
    robo.FV = [var('FV{0}'.format(i)) for i in robo.num]
    robo.MS = [Matrix(var('MX{0}, MY{0}, MZ{0}'.format(i))) for i in robo.num]
    robo.M = [var('M{0}'.format(i)) for i in robo.num]
    robo.GAM = [var('GAM{0}'.format(i)) for i in robo.num]
    inertia_matrix_terms = ("XX{0}, XY{0}, XZ{0}, ") + \
        ("XY{0}, YY{0}, YZ{0}, ") + \
        ("XZ{0}, YZ{0}, ZZ{0}")
    robo.J = [
        Matrix(3, 3, var(inertia_matrix_terms.format(i))) \
        for i in robo.num
    ]
    robo.G = Matrix([0, 0, -var('G3')])
    robo.w0 = zeros(3, 1)
    robo.wdot0 = zeros(3, 1)
    robo.v0 = zeros(3, 1)
    robo.vdot0 = zeros(3, 1)
    robo.q = [0, var('q1'), var('q2')]
    robo.qdot = [0, var('QP1'), var('QP2')]
    robo.qddot = [0, var('QDP1'), var('QDP2')]
    return robo


def sr400():
    #TODO: bring it to the new notation with 0-frame
    """Generate Robot instance of SR400"""
    robo = Robot('SR400', 8, 9, 10, False)
    robo.ant = [-1, 0, 1, 2, 3, 4, 5, 1, 7, 8, 3]
    robo.sigma = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2]
    robo.mu = [0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
    robo.alpha = [0, 0, -pi/2, 0, -pi/2, pi/2, -pi/2, -pi/2, 0, 0, 0]
    d_var = var('D:9')
    robo.d = [0, 0, d_var[2], d_var[3], d_var[4], 0, 0,
              d_var[2], d_var[8], d_var[3], -d_var[8]]
    robo.theta = [0] + list(var('th1:10')) + [0]
    robo.r = [0, 0, 0, 0, var('RL4'), 0, 0, 0, 0, 0, 0]
    robo.b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    robo.gamma = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pi/2]
    robo.structure = tools.CLOSED_LOOP
    return robo

def SimpleLoop():
    """Generate Robot instance of RRRR loop"""
    robo = Robot('SR400', 3, 4, 5, False)
    robo.ant = [-1, 0, 1, 0, 3, 2]
    robo.sigma = [0, 0, 0, 0, 0, 2]
    robo.mu = [0, 1, 0, 0, 0, 0]
    robo.alpha = [0, 0, 0, 0, 0, 0]
    d_var = var('D:9')
    robo.d = [0, 0, d_var[2], d_var[3], d_var[3], d_var[2]]
    robo.theta = [0] + list(var('th1:5')) + [0]
    robo.r = [0, 0, 0, 0, 0, 0]
    robo.b = [0, 0, 0, 0, 0, 0]
    robo.gamma = [0, 0, 0, 0, 0, 0]
    robo.structure = tools.CLOSED_LOOP
    return robo

def SimpleLoop2():
    """Generate Robot instance of PRRR loop"""
    robo = Robot('LOL', 3, 4, 5, False)
    robo.ant = [-1, 0, 0, 2, 3, 1]
    robo.sigma = [0, 1, 0, 0, 0, 2]
    robo.mu = [0, 1, 0, 0, 0, 0]
    robo.alpha = [0, 0, pi/2, 0, 0, pi/2]
    d_var = var('D:9')
    robo.d = [0, 0, 0, d_var[3], d_var[3],0]
    robo.theta = [0, 0] + list(var('th2:5')) + [0]
    robo.r = [0, var('r'), 0, 0, 0, 0]
    robo.b = [0, 0, 0, 0, 0, 0]
    robo.gamma = [0, 0, 0, 0, 0, 0]
    robo.structure = tools.CLOSED_LOOP
    return robo

def rx90():
    """Generate Robot instance of RX90"""
    robo = Robot('RX90', 6, 6, 6, False)
    # table of geometric parameters RX90
    robo.sigma = [2, 0, 0, 0, 0, 0, 0]
    robo.alpha = [0, 0, pi/2, 0, -pi/2, pi/2, -pi/2]
    robo.d = [0, 0, 0, var('D3'), 0, 0, 0]
    robo.theta = [0] + list(var('th1:7'))
    robo.r = [0, 0, 0, 0, var('RL4'), 0, 0]
    robo.b = [0, 0, 0, 0, 0, 0, 0]
    robo.gamma = [0, 0, 0, 0, 0, 0, 0]
    robo.mu = [0, 1, 1, 1, 1, 1, 1]
    robo.structure = tools.SIMPLE
    robo.w0 = zeros(3, 1)
    robo.wdot0 = zeros(3, 1)
    robo.v0 = zeros(3, 1)
    robo.vdot0 = zeros(3, 1)
    num = range(0, 7)
    robo.qdot = [var('QP{0}'.format(i)) for i in num]
    robo.qddot = [var('QDP{0}'.format(i)) for i in num]
    robo.Nex= [zeros(3, 1) for i in num]
    robo.Nex[-1] = Matrix(var('CX{0}, CY{0}, CZ{0}'.format(robo.NJ)))
    robo.Fex = [zeros(3, 1) for i in num]
    robo.Fex[-1] = Matrix(var('FX{0}, FY{0}, FZ{0}'.format(robo.NJ)))
    robo.FS = [var('FS{0}'.format(i)) for i in num]
    robo.IA = [var('IA{0}'.format(i)) for i in num]
    robo.FV = [var('FV{0}'.format(i)) for i in num]
    robo.MS = [Matrix(var('MX{0}, MY{0}, MZ{0}'.format(i))) for i in num]
    robo.M = [var('M{0}'.format(i)) for i in num]
    robo.GAM = [var('GAM{0}'.format(i)) for i in num]
    inertia_matrix_terms = ("XX{0}, XY{0}, XZ{0}, ") + \
        ("XY{0}, YY{0}, YZ{0}, ") + \
        ("XZ{0}, YZ{0}, ZZ{0}")
    robo.J = [
        Matrix(3, 3, var(inertia_matrix_terms.format(i))) \
        for i in num
    ]
    robo.G = Matrix([0, 0, var('G3')])
    return robo



