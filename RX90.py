# -*- coding: utf-8 -*-
"""
Created on Sun Jun 02 09:42:00 2013

@author: Bogdqn
"""
from sympy import pi,var,zeros,Matrix
D3 = var('D3')
RL4 = var('RL4')
#table of geometric parameters RX90
NJ = 6
NL = 6
num = range(1,NJ+1)
ant = (-1,0,1,2,3,4)
sigm = (0,0,0,0,0,0)
alpha = (0,pi/2,0,-pi/2,pi/2,-pi/2)
d = (0,0,D3,0,0,0)
theta = list(var('th1:7'))
r = (0,0,0,RL4,0,0)
b = (0,0,0,0,0,0)
gamma = (0,0,0,0,0,0)

W0 = zeros(3,1)
WP0 = zeros(3,1)
V0 = zeros(3,1)
VP0 = zeros(3,1)

q = [var('Q{0}'.format(i)) for i in num]
qp =  [var('QP{0}'.format(i)) for i in num]
qdp = [var('QDP{0}'.format(i)) for i in num]

nej = [zeros(3,1) for i in num]
nej[-1] = Matrix(var('CX{0},CY{0},CZ{0}'.format(num[-1])))
fej = [zeros(3,1) for i in num]
fej[-1] = Matrix(var('FX{0},FY{0},FZ{0}'.format(num[-1])))
FS = var('FS0:{0}'.format(NL))
IA = var('IA0:{0}'.format(NL))
FV = var('FV0:{0}'.format(NL))
MS = [Matrix(var('MX{0},MY{0},MZ{0}'.format(i))) for i in num]
M = [var('M{0}'.format(i)) for i in num]
GAM = [0 for i in num]
J = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                            'XY{0},YY{0},YZ{0},'
                            'XZ{0},YZ{0},ZZ{0}').format(i))) for i in num]
G = Matrix([0,0,var('G3')])