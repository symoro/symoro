# -*- coding: utf-8 -*-
"""
Created on Sun Jun 02 09:42:00 2013

@author: Bogdqn
"""
from sympy import pi,var
D3 = var('D3')
RL4 = var('RL4')
#table of geometric parameters RX90
ant = (-1,0,1,2,3,4,5)
sigm = (2,0,0,0,0,0,0)
alpha = (0,0,pi/2,0,-pi/2,pi/2,-pi/2)
d = (0,0,0,D3,0,0,0)
theta = var('th:7')
r = (0,0,0,0,RL4,0,0)
b = (0,0,0,0,0,0,0)
gamma = (0,0,0,0,0,0,0)
NJ = len(sigm)
