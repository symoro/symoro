# -*- coding: utf-8 -*-
"""
Created on Sun Jun 02 10:24:44 2013

@author: Bogdqn
"""

from sympy import pi,var,zeros, Matrix
#table of geometric parameters RX90
ant = (-1,0)
sigm = (1,0)
alpha = (pi/2,pi/2)
d = (0,0)
theta = (pi/2,var('th2'))
r = (var('x'),0)
b = (0,0)
gamma = (0,0)
num = range(1,3)
NJ = 2
NL = 2
n_ext = [zeros(3,1) for i in num]
f_ext = [zeros(3,1) for i in num]
FS = [0 for i in num]
IA = [0 for i in num]
FV = [0 for i in num]
MS = [zeros(3,1) for i in num]
MS[1][0] = var('m2l')
M = [var('m{0}'.format(i)) for i in num]
GAM = [var('u'),0]
J = [zeros(3) for i in num]
J[1][2,2] = var('m2ll')
G = Matrix([0,0,-var('g')])
W0 = zeros(3,1)
WP0 = zeros(3,1)
V0 = zeros(3,1)
VP0 = zeros(3,1)
q = var('x,th')
qp =  var('xd,thd')
qdp = var('xdd,thdd')