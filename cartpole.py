# -*- coding: utf-8 -*-
"""
Created on Sun Jun 02 10:24:44 2013

@author: Bogdqn
"""

from sympy import pi,var
#table of geometric parameters RX90
ant = (-1,0,1)
sigm = (2,1,0)
alpha = (0,pi/2,pi/2)
d = (0,0,0)
theta = (0,pi/2,var('th2'))
r = (0,var('r1'),0)
b = (0,0,0)
gamma = (0,0,0)
NJ = len(sigm)