# -*- coding: utf-8 -*-
"""
Created on Mon May 13 21:35:46 2013

@author: Bogdqn
"""


from sympy import *

def rot_trans(axis='z', th=0, p=0):
    if axis == 'x':
        return  Matrix([[1, 0, 0, p],
                        [0, cos(th),-sin(th), 0],
                        [0, sin(th), cos(th), 0],
                        [0, 0, 0, 1]])
    elif axis == 'y':
        return  Matrix([[cos(th), 0, sin(th), 0],
                        [0, 1, 0, p],
                        [-sin(th),0, cos(th), 0],
                        [0, 0, 0, 1]])
    else:
        return  Matrix([[cos(th),-sin(th), 0,0],
                        [sin(th),cos(th), 0,0],
                        [0,0, 1,p],
                        [0,0,0,1]])

def transform(theta, r, alpha, d, gamma=0, b=0):
    return rot_trans('z',gamma,b)*rot_trans('x',alpha,d)*rot_trans('z',theta,r)

def transform_inv(theta, r, alpha, d, gamma=0, b=0):
    return rot_trans('z',-theta,-r)*rot_trans('x',-alpha,-d)*rot_trans('z',-gamma,-b)

def sym_replace(old_sym, sym_dict, name):
    new_sym = symbols(name)
    sym_dict[new_sym] = old_sym
    return new_sym

def mat_trigsimp(M):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            M[i1,i2] = trigsimp(M[i1,i2])

def trig_replace(M,sym_dict,angle,number):
    sin_sym,cos_sym = symbols('S{0},C{0}'.format(number))
    #save sin and cos to the dictionaty
    sym_dict[cos_sym] = cos(angle)
    sym_dict[sin_sym] = sin(angle)
    return M.subs({sym_dict[sin_sym]:sin_sym, sym_dict[cos_sym]:cos_sym})

def mat_trig_replace(M,sym_dict,index_list,theta):
    number = ''
    angle = sympify(0)
    for index in index_list:
        M = trig_replace(M,sym_dict,theta[index],index+1)
        number = '{0}{1}'.format(index+1,number)
        angle += theta[index]
    return trig_replace(M,sym_dict,angle,number)

def mat_sym_replace(M, sym_dict, name='M'):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            if M[i1,i2].count_ops() != 0 and (-M[i1,i2]).count_ops() != 0:
                M[i1,i2] = sym_replace(M[i1,i2],sym_dict, '{0}{1}{2}'.format(name,i1+1,i2+1))
                print M[i1,i2], '=', sym_dict[M[i1,i2]]

#needed variables
D3 = symbols('D3')
RL4 = symbols('RL4')
#table of geometric parameters RX90
ant = (0,1,2,3,4,5)
sigm = (0,0,0,0,0,0)
alpha = (0,pi/2,0,-pi/2,pi/2,-pi/2)
d = (0,0,D3,0,0,0)
theta = var('th1,th2,th3,th4,th5,th6')
r = (0,0,0,RL4,0,0)
b = (0,0,0,0,0,0)
gamma = (0,0,0,0,0,0)
#identity matrix as initial transformation
T_prev = eye(4)
#dictionary for symbol substitution
sym_dict = {}
trig_dict = {}
#computing DGM
i = 5
i_term = -1
j_last = i+1
T_list = {}
while i != i_term:
    index = 0
    index_list = []
    T = eye(4)
    while True:
        T = transform(theta[i],r[i],alpha[i],d[i],gamma[i],b[i])*T
        mat_trigsimp(T)
        index = index*10 + i+1
        index_list.append(i)
        i_prev = i
        i = ant[i]-1
        if gamma[i_prev] != 0 or alpha[i_prev] != 0 or i == i_term:
            break

    T = mat_trig_replace(T,trig_dict,index_list,theta)
    T_next = T*T_prev
    #make substitution with saving the symbols to the dictionary
    mat_sym_replace(T_next,sym_dict,'T{0}{1}'.format(i+1,j_last))
    T_list['T{0}{1}'.format(i+1,j_last)] = T_next
    T_prev = T_next

#making the substitution
print 'unfolded DGM:'
for i2 in range(4):
    for i1 in range(3):
        print 'T0{2}{0}{1}'.format(i1+1,i2+1,j_last), '=', T_next[i1,i2].subs(sym_dict)
