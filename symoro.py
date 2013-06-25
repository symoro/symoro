# -*- coding: utf-8 -*-
"""
Created on Wed May 29 19:57:45 2013

@author: Bogdqn
"""
from sympy import Matrix, eye, zeros, sin, cos, Integer, var, trigsimp, sympify

axis_dict = {'x':0,'y':1,'z':2}
           
def get_r(T):
    return T[:3,:3]

def get_p(T):
    return T[:3,3]

#????????????????????????????????
def screw_transform(T):
    R, P = get_rp(T)
    return screw_transf_rp(R,P)

#??????????????????????????????????
def screw_transf_rp(R,P):
    return Matrix([R.row_join(-R*hat(P)),zeros(3).row_join(R)])

def hat(v):
    return Matrix([[0,-v[2],v[1]],
                   [v[2],0,-v[0]],
                   [-v[1],v[0],0]])

def rot(axis = 'z', th = 0):
    if axis == 'x':
        return  Matrix([[1, 0, 0],
                        [0, cos(th),-sin(th)],
                        [0, sin(th), cos(th)]])
    elif axis == 'y':
        return  Matrix([[cos(th), 0, sin(th)],
                        [0, 1, 0],
                        [-sin(th),0, cos(th)]])
    else:
        return  Matrix([[cos(th),-sin(th), 0],
                        [sin(th),cos(th), 0],
                        [0,0, 1]])

def trans_vect(axis = 'z', p = 0):
    v = zeros(3,1)
    v[axis_dict[axis]] = p
    return v

def trans(axis = 'z', p = 0):
    return Matrix([eye(3).row_join(trans_vect(axis,p)),
                       [0,0,0,1]])

def rot_trans(axis='z', th=0, p=0):
    return Matrix([rot(axis,th).row_join(trans_vect(axis,p)),
                       [0,0,0,1]])

def transform(theta, r, alpha, d, gamma=0, b=0, invert = False):
    if not invert:
        return rot_trans('z',gamma,b)*rot_trans('x',alpha,d)*rot_trans('z',theta,r)
    else:
        return rot_trans('z',-theta,-r)*rot_trans('x',-alpha,-d)*rot_trans('z',-gamma,-b)

def mat_trigsimp(M):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            M[i1,i2] = trigsimp(M[i1,i2])

def trig_replace(M,sym_dict,angle,number,disp = True):
#    print number
    sin_sym,cos_sym = var('S{0},C{0}'.format(number))
    sym_dict[cos_sym] = cos(angle)
    sym_dict[sin_sym] = sin(angle)
    if disp:
        print cos_sym, '=', sym_dict[cos_sym]
        print sin_sym, '=', sym_dict[sin_sym]
    return M.subs({sym_dict[sin_sym]:sin_sym, sym_dict[cos_sym]:cos_sym})


def mat_trig_replace(M,sym_dict,index_list,theta,alpha,gamma,disp = True):
    number = ''
    angle = Integer(0)
    for j in index_list:
        index = j + 1
        if not sympify(theta[j]).is_constant():
            M = trig_replace(M,sym_dict,theta[j],index,disp)
        if not sympify(alpha[j]).is_constant():
            M = trig_replace(M,sym_dict,alpha[j],index,disp)
        if not sympify(gamma[j]).is_constant():
            M = trig_replace(M,sym_dict,gamma[j],index,disp)
            if alpha[j] == 0:
                number = 'G' + str(index) + number
                angle += gamma[j]
        number = str(index)+number
        angle += theta[j]
        if angle != theta[j] and not sympify(angle).is_constant():
            M = trig_replace(M,sym_dict,angle,number)
    return M

def sym_replace(old_sym,sym_dict,name,index='',disp=True,forced=False):
    if old_sym.count_ops() != 0 and (-old_sym).count_ops() != 0 or forced:
        if index != None:
            name_index = name+str(index)
        new_sym = var(name_index)
        sym_dict[new_sym] = old_sym
        if disp:
            print new_sym, '=', old_sym
        return new_sym
    else:
        return old_sym

def mat_sym_replace(M,sym_dict,name='MATRIX',index = '',disp = True):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            if M[i1,i2].count_ops() != 0 and (-M[i1,i2]).count_ops() != 0:
                if M.shape[1] > 1:
                    name_index = name+str(i1+1)+str(i2+1)
                else:
                    name_index = name+str(i1+1)
                M[i1,i2] = sym_replace(M[i1,i2],sym_dict, name_index,index,disp = disp)
    return M
    
def matsymm_sym_replace(M, sym_dict, name='MATRIX',index = '',disp = True):
    for i2 in range(M.shape[1]):
        for i1 in range(i2,M.shape[0]):
            name_index = name+str(i1+1)+str(i2+1)
            M[i1,i2] = sym_replace(M[i1,i2],sym_dict, name_index,index,disp = disp)
            M[i2,i1] = M[i1,i2]
    return M

def unfold(expr,symb_dict):
    while any([symb_dict.has_key(a) for a in expr.atoms()]):
        expr = expr.subs(symb_dict)
    return expr

def mat_unfold(M,symb_dict):
    for i1 in range(M.shape[0]):
        for i2 in range(M.shape[1]):
            M[i1,i2] = unfold(M[i1,i2],symb_dict)

#P = Matrix(var('PX,PY,PZ'))
#ang = var('theta,R,alpha,D,gamma,B')
#T = transform(ang[0],ang[1],ang[2],ang[3],ang[4],ang[5],True)
#print get_r(T)*P + get_p(T)
#print transform(ang[0],ang[1],ang[2],ang[3],ang[4],ang[5])#*rot_trans(axis = 'x',th = var('a2'))