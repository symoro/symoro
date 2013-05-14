# -*- coding: utf-8 -*-
"""
Created on Mon May 13 21:35:46 2013

@author: Bogdqn
"""

from sympy import *

#needed variables
D3 = symbols('D3')
RL4 = symbols('RL4')
#table of geometric parameters RX90
num = range(1,7)
sigm = (0,0,0,0,0,0)
alfa = (0,pi/2,0,-pi/2,pi/2,-pi/2)
d = (0,0,D3,0,0,0)
theta = var('th:1:7')
r = (0,0,0,RL4,0,0)

#identity matrix as initial transformation
Tprev = Matrix(4,4,lambda i,j: int(i==j))
#dictionary for symbol substitution
sym_list = {}

#computing DGM
j_last = len(num)-1
for i in range(j_last,-1,-1):
    j = num[i]
    th  = theta[i]
    al = alfa[i]

    cos_sym,sin_sym = symbols('S{0},C{0}'.format(j))
    #save sin and cos to the dictionaty
    sym_list[cos_sym] = cos(th)
    sym_list[sin_sym] = sin(th)
    T = Matrix([[cos_sym,-sin_sym, 0, d[i]],
                [cos(al)*sin_sym,cos(al)*cos_sym, -sin(al),-r[i]*sin(al)],
                [sin(al)*sin_sym,sin(al)*cos_sym, cos(al),r[i]*cos(al)],
                [0,0,0,1]])
    Tnext = T*Tprev

    #make substitution with saving the symbols to the dictionary
    for i1 in range(4):
        for i2 in range(4):
            if Tnext[i1,i2] != 0:
                new_sym = symbols('T{0}{1}{2}{3}'.format(j,j_last,i1+1,i2+1))
                sym_list[new_sym] = Tnext[i1,i2]
                Tnext[i1,i2] = new_sym
                print new_sym, '=', sym_list[new_sym]

    Tprev = Tnext

#making the substitution
print Tprev.subs(sym_list)
