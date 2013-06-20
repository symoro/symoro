# -*- coding: utf-8 -*-
"""
Created on Sun Jun 02 10:24:44 2013

@author: Bogdqn
"""

from sympy import pi,var,zeros,Matrix
#table of geometric parameters RX90
ant = (-1,0)
sigm = (0,0)
alpha = (0,0)
d = var('L1:3')
theta = var('th1:3')
r = (0,0)
b = (0,0)
gamma = (0,0)
num = range(1,3)
NJ = 2
NL = 2
G = Matrix([var('g'),0,0])
nej = [zeros(3,1) for i in num]
fej = [zeros(3,1) for i in num]
FS = [0 for i in num]
IA = [0 for i in num]
FV = [0 for i in num]
MS = [zeros(3,1) for i in num]
MS[0][0] = var('m1l')
MS[1][0] = var('m2l')
M = [var('m{0}'.format(i)) for i in num]
GAM = [0,var('u')]
J = [zeros(3) for i in num]
J[0][2,2] = var('I1')
J[1][2,2] = var('I2')
W0 = zeros(3,1)
WP0 = zeros(3,1)
V0 = zeros(3,1)
VP0 = zeros(3,1)
q = theta
qp =  var('thd1:3')
qdp = var('thdd1:3')