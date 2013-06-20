# -*- coding: utf-8 -*-
"""
Created on Sat Jun 01 23:13:36 2013

@author: Bogdqn
"""

from sympy import sign, zeros
from symoro import *
from cartpole import *

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
    wp[j] =  mat_sym_replace(jAant*wpant+(1-sigm[j])*(qdpj+hat(wi[j])*qpj),
                            sydi,'WP',index)                   
    U[j] = mat_sym_replace(hat(w[j])*hat(w[j])+hat(wp[j]),sydi,'U',index)
    vsp[j] = mat_sym_replace(vpant + Uant*antPj[j],sydi,'VSP',index)
    vp[j] = mat_sym_replace(jAant*vsp[j] + sigm[j]*(qdpj+2*hat(wi[j])*qpj),
                                sydi,'VP',index)                         
    F[j] = mat_sym_replace(M[j]*vp[j] + U[j]*MS[j],sydi,'F',index)
    Psi[j] = mat_sym_replace(J[j]*w[j],sydi,'PSI',index)
    N[j] = mat_sym_replace(J[j]*wp[j] + hat(w[j])*Psi[j],
                        sydi,'No',index)
for j in reversed(range(NL)):
    index = num[j]
    fj[j] = mat_sym_replace(F[j] + fej[j],sydi,'E',index)
    nj[j] = mat_sym_replace(N[j] + nej[j] + hat(MS[j])*vp[j],sydi,'N',index)
    i = ant[j]
    fi[j] = mat_sym_replace(antAj[j]*fj[j],sydi,'FDI',index)
    if i != -1:
        fej[i] += fi[j]
        nej[i] += antAj[j]*nj[j] + hat(antPj[j])*fi[j]
for j in range(NL):
    index = num[j]
    if sigm[j] != 2:
        GAM[j] = sym_replace((sigm[j]*fj[j]+(1-sigm[j])*nj[j])[2]+FS[j]*sign(qp[j])+
                                FV[j]*qp[j]+IA[j]*qdp[j],sydi,'GAM',index,forced = True)

#for j in range(NL):
#    print  'GAM{0}'.format(num[j]), '=', unfold(GAM[j],sydi)

#    DV11 = sym_replace(-w[j][0]*w[j][0],sydi,'DV11',index)
#    DV22 = sym_replace(-w[j][1]*w[j][1],sydi,'DV22',index)
#    DV33 = sym_replace(-w[j][2]*w[j][2],sydi,'DV33',index)
#    DV12 = sym_replace(w[j][0]*w[j][1],sydi,'DV12',index)
#    DV13 = sym_replace(w[j][0]*w[j][2],sydi,'DV13',index)
#    DV23 = sym_replace(w[j][1]*w[j][2],sydi,'DV23',index)
#    Matrix([[DV33+DV22,DV12,DV13],
#                                  [DV12,DV11+DV33,DV23],
#                                    [DV13,DV23,DV11+DV22]])