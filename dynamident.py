# -*- coding: utf-8 -*-
"""
Created on Sat Jun 08 07:35:39 2013

@author: Bogdqn
"""
# -*- coding: utf-8 -*-
"""
Created on Sat Jun 01 23:13:36 2013

@author: Bogdqn
"""

from sympy import sign, zeros
from symoro import *
from RX90 import *

w = [zeros(3,1) for i in num]
wi = [zeros(3,1) for i in num]
wp = [zeros(3,1) for i in num]
v = [zeros(3,1) for i in num]
vp = [zeros(3,1) for i in num]
vsp = [zeros(3,1) for i in num]
phij = [zeros(6,1) for i in num]
lamdaj = [0 for i in num]
effappj = [zeros(3,1) for i in num]
antAj = [zeros(3,3) for i in num]
antPj = [zeros(3,1) for i in num]
U = [zeros(3,3) for i in num]
F = [zeros(3,1) for i in num]
Psi = [zeros(3,1) for i in num]
N = [zeros(3,1) for i in num]
fj = [zeros(3,1) for i in num]
fi = [zeros(3,1) for i in num]
nj = [zeros(3,1) for i in num]
GAM = [zeros(NL,13) for i in num]
sydi = {}
tridi = {}
for j in range(NL):
    index = num[j]
    antTj = mat_trig_replace(transform(theta[j],r[j],alpha[j],d[j]),tridi,[j],theta,alpha,gamma)
    antAj[j] = mat_sym_replace(get_r(antTj),sydi,'A',index)
    antPj[j] = mat_sym_replace(get_p(antTj),sydi,'L',index)
for j in range(NL):
    index = num[j]
    jAant = antAj[j].T
    qpj = Matrix([0,0,qp[j]])
    qdpj = Matrix([0,0,qdp[j]])
    if j != 0:
        Uant = U[ant[j]]
        want = w[ant[j]]
        wpant = wp[ant[j]]
        vpant = vp[ant[j]]
    else:
        Uant = zeros(3,3)
        want = W0
        wpant = WP0
        vpant = VP0-G
    wi[j] = mat_sym_replace(jAant*want,sydi,'WI',index)
    w[j] = mat_sym_replace(wi[j]+(1-sigm[j])*qpj,sydi,'W',index)
    wp[j] = mat_sym_replace(jAant*wpant+(1-sigm[j])*(qdpj+hat(wi[j])*qpj),
                            sydi,'WP',index)                   
    U[j] = mat_sym_replace(hat(w[j])*hat(w[j])+hat(wp[j]),sydi,'U',index)
    vsp[j] = mat_sym_replace(vpant + Uant*antPj[j],sydi,'VSP',index)
    vp[j] = mat_sym_replace(jAant*vsp[j] + sigm[j]*(qdpj+2*hat(wi[j])*qpj),
                                sydi,'VP',index)                         
 
def get_dyn_param(P):
    J = Matrix([[P[0],P[1],P[2]],
                [P[1],P[3],P[4]],
                [P[2],P[4],P[5]]])
    MS = Matrix(3,1,P[6:9])
    M = P[9]
    return J,MS,M
    
for k in range(NL):
    param_vec = [J[k][0],J[k][1],J[k][2],J[k][4],J[k][5],J[k][8],
                 MS[k][0],MS[k][1],MS[k][2],M[k]]
    for i in range(10):
        if param_vec[i] == Integer(0):
            continue
        index = str(num[k])+str(param_vec[i])
        mask = zeros(10,1)
        mask[i] = 1
        J_tmp,MS_tmp,M_tmp = get_dyn_param(mask)
        F[j] = mat_sym_replace(M_tmp*vp[j] + U[j]*MS_tmp,sydi,'F',index)
        Psi[j] = mat_sym_replace(J_tmp*w[j],sydi,'PSI',index)
        N[j] = mat_sym_replace(J_tmp*wp[j] + hat(w[j])*Psi[j] + hat(MS_tmp)*vp[j],
                            sydi,'No',index)
                        
        for j in reversed(range(k+1)):
            index = str(num[j])+str(param_vec[i])
            fj[j] = mat_sym_replace(F[j] + fej[j],sydi,'E',index)
            nj[j] = mat_sym_replace(N[j] + nej[j],sydi,'N',index)
            aj = ant[j]
            if aj != -1:
                fej[aj] = antAj[j]*fj[j]
                nej[aj] = antAj[j]*nj[j] + hat(antPj[j])*fi[j]
        for j in range(k+1):
            index = str(num[j])+str(param_vec[i])
            if sigm[j] != 2:
                GAM[j][k,i] = sym_replace((sigm[j]*fj[j]+(1-sigm[j])*nj[j])[2],sydi,'DG',index,forced = True)
        
    if IA[k] != Integer(0) and qdp[k] != Integer(0):
        index = str(num[k])+str(IA[k])
        GAM[k][k,10] = sym_replace(qdp[k],sydi,'DG',index,forced = True)
        
    if FS[k] != Integer(0) and qp[k] != Integer(0):
        index = str(num[k])+str(IA[k])
        GAM[k][k,10] = sym_replace(sign(qp[k]),sydi,'DG',index,forced = True)
        
    if FV[k] != Integer(0) and qp[k] != Integer(0):
        index = str(num[k])+str(IA[k])
        GAM[k][k,10] = sym_replace(qp[k],sydi,'DG',index,forced = True)