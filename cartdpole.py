# -*- coding: utf-8 -*-
"""
Created on Sun Jun 02 10:24:44 2013

@author: Bogdqn
"""

from sympy import pi,var,zeros, Matrix
#table of geometric parameters RX90
ant = (-1,0,0)
sigm = (1,0,0)
alpha = (pi/2,pi/2,pi/2)
d = (0,0,0)
theta = (pi/2,var('th[1]'),var('th[2]'))
r = (var('th[0]'),0,0)
b = (0,0,1)
gamma = (0,0,0)
num = range(1,4)
NJ = 3
NL = 3
nej = [zeros(3,1) for i in num]
fej = [zeros(3,1) for i in num]
FS = [0 for i in num]
IA = [0 for i in num]
FV = [0 for i in num]
MS = [zeros(3,1) for i in num]
MS[1][0] = var('m2l')
MS[2][0] = var('m2l')
M = [var('m{0}'.format(i)) for i in num]
GAM = [var('u[0]'),var('u[1]'),0]
J = [zeros(3) for i in num]
J[1][2,2] = var('m2ll')
J[2][2,2] = var('m2ll')
G = Matrix([0,0,-var('g')])
W0 = zeros(3,1)
WP0 = zeros(3,1)
V0 = zeros(3,1)
VP0 = zeros(3,1)
q = [r[0],theta[1],theta[2]]
qp =  var('thd[:3]')
qdp = var('thdd[:3]')