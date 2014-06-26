# -*- coding: utf-8 -*-


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
    # TODO: bring it to the new notation with 0-frame
    robot = Robot()
    robot.name = 'CartPole'
    robot.ant = (-1, 0)
    robot.sigma = (1, 0)
    robot.alpha = (pi/2, pi/2)
    robot.d = (0, 0)
    robot.theta = (pi/2, var('Th2'))
    robot.r = (var('R1'), 0)
    robot.b = (0, 0)
    robot.gamma = (0, 0)
    robot.num = range(1, 3)
    robot.NJ = 2
    robot.NL = 2
    robot.NF = 2
    robot.Nex = [zeros((3, 1))] * robot.num
    robot.Fex = [zeros((3, 1))] * robot.num
    robot.FS = [0] * robot.num
    robot.IA = [0] * robot.num
    robot.FV = [var('FV{0}'.format(i)) for i in robot.num]
    robot.MS = [zeros((3, 1))] * robot.num
    robot.MS[1][0] = var('MX2')
    robot.M = [var('M{0}'.format(i)) for i in robot.num]
    robot.GAM = [var('GAM{0}'.format(i)) for i in robot.num]
    robot.J = [zeros(3)] * robot.num
    robot.J[1][2, 2] = var('ZZ2')
    robot.G = Matrix([0, 0, -var('G3')])
    robot.w0 = zeros((3, 1))
    robot.wdot0 = zeros((3, 1))
    robot.v0 = zeros((3, 1))
    robot.vdot0 = zeros((3, 1))
    robot.q = var('R1, Th2')
    robot.qdot = var('R1d, Th2d')
    robot.qddot = var('R1dd, Th2dd')
    robot.num.append(0)
    return robot


def planar2r():
    """Generate Robot instance of 2R Planar robot"""
    robot = Robot('Planar2R', 2, 2, 3, False)
    robot.structure = tools.SIMPLE
    robot.sigma = [2, 0, 0, 2]
    robot.mu = [0, 1, 1, 0]
    robot.gamma = [0, 0, 0, 0]
    robot.b = [0, 0, 0, 0]
    robot.alpha = [0, 0, 0, 0]
    robot.d = [0, 0, var('L1'), var('L2')]
    robot.theta = [0, var('th1'), var('th2'), 0]
    robot.r = [0, 0, 0, 0]
    return robot


def sr400():
    #TODO: bring it to the new notation with 0-frame
    """Generate Robot instance of SR400"""
    robot = Robot('SR400', 8, 9, 10, False)
    robot.ant = [-1, 0, 1, 2, 3, 4, 5, 1, 7, 8, 3]
    robot.sigma = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 2]
    robot.mu = [0, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0]
    robot.alpha = [0, 0, -pi/2, 0, -pi/2, pi/2, -pi/2, -pi/2, 0, 0, 0]
    d_var = var('D:9')
    robot.d = [0, 0, d_var[2], d_var[3], d_var[4], 0, 0,
              d_var[2], d_var[8], d_var[3], -d_var[8]]
    robot.theta = [0] + list(var('th1:10')) + [0]
    robot.r = [0, 0, 0, 0, var('RL4'), 0, 0, 0, 0, 0, 0]
    robot.b = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    robot.gamma = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, pi/2]
    robot.structure = tools.CLOSED_LOOP
    return robot


def rx90():
    """Generate Robot instance of RX90"""
    robot = Robot('RX90', 6, 6, 6, False)
    # table of geometric parameters RX90
    robot.sigma = [2, 0, 0, 0, 0, 0, 0, 0]
    robot.alpha = [0, 0, pi/2, 0, -pi/2, pi/2, -pi/2]
    robot.d = [0, 0, 0, var('D3'), 0, 0, 0]
    robot.theta = [0] + list(var('th1:7'))
    robot.r = [0, 0, 0, 0, var('RL4'), 0, 0]
    robot.b = [0, 0, 0, 0, 0, 0, 0]
    robot.gamma = [0, 0, 0, 0, 0, 0, 0]
    robot.mu = [0, 1, 1, 1, 1, 1, 1]
    robot.structure = tools.SIMPLE
    # robot.w0 = zeros((3, 1))
    # robot.wdot0 = zeros((3, 1))
    # robot.v0 = zeros((3, 1))
    # robot.vdot0 = zeros((3, 1))
    # robot.qdot = [var('QP{0}'.format(i)) for i in num]
    # robot.qddot = [var('QDP{0}'.format(i)) for i in num]
    # robot.Nex= [zeros((3, 1)) for i in num]
    # robot.Nex[-1] = Matrix(var('CX{0}, CY{0}, CZ{0}'.format(robot.NJ)))
    # robot.Fex = [zeros((3, 1)) for i in num]
    # robot.Fex[-1] = Matrix(var('FX{0}, FY{0}, FZ{0}'.format(robot.NJ)))
    # robot.FS = [var('FS{0}'.format(i)) for i in num]
    # robot.IA = [var('IA{0}'.format(i)) for i in num]
    # robot.FV = [var('FV{0}'.format(i)) for i in num]
    # robot.MS = [Matrix(var('MX{0}, MY{0}, MZ{0}'.format(i))) for i in num]
    # robot.M = [var('M{0}'.format(i)) for i in num]
    # robot.GAM = [var('GAM{0}'.format(i)) for i in num]
    # robot.J = [Matrix(3, 3, var(('XX{0}, XY{0}, XZ{0}, '
    #                     'XY{0}, YY{0}, YZ{0}, '
    #                     'XZ{0}, YZ{0}, ZZ{0}').format(i))) for i in num]
    # robot.G = Matrix([0, 0, var('G3')])
    # robot.num.append(0)
    return robot
