# -*- coding: utf-8 -*-
"""
Created on Sat Jun 01 23:13:36 2013

@author: Bogdqn
"""

from sympy import pi,sign,sympify,simplify
from symoro import *
from cartpole import *

q = var('Q0:{0}'.format(NJ))
qp = var('QP0:{0}'.format(NJ))
qdp = var('QDP0:{0}'.format(NJ))
w = [zeros(3,1) for i in range(NJ)]
wi = [zeros(3,1) for i in range(NJ)]
wp = [zeros(3,1) for i in range(NJ)]
v = [zeros(3,1) for i in range(NJ)]
vp = [zeros(3,1) for i in range(NJ)]
vsp = [zeros(3,1) for i in range(NJ)]
phij = [zeros(6,1) for i in range(NJ)]
lamdaj = [0 for i in range(NJ)]
effappj = [zeros(3,1) for i in range(NJ)]
antAj = [zeros(3,3) for i in range(NJ)]
antPj = [zeros(3,1) for i in range(NJ)]
U = [zeros(3,3) for i in range(NJ)]
F = [zeros(3,1) for i in range(NJ)]
Psi = [zeros(3,1) for i in range(NJ)]
N = [zeros(3,1) for i in range(NJ)]
fj = [zeros(3,1) for i in range(NJ)]
fi = [zeros(3,1) for i in range(NJ)]
nj = [zeros(3,1) for i in range(NJ)]
nej = [zeros(3,1) for i in range(NJ)]#[Matrix(var('CX{0},CY{0},CZ{0}'.format(i))) for i in range(NJ)]
fej = [zeros(3,1) for i in range(NJ)]#[Matrix(var('FX{0},FY{0},FZ{0}'.format(i))) for i in range(NJ)]
FS = var('FS0:{0}'.format(NJ))
IA = var('IA0:{0}'.format(NJ))
FV = var('FV0:{0}'.format(NJ))
MS = [zeros(3,1) for i in range(NJ)]#[Matrix(var('MX{0},MY{0},MZ{0}'.format(i))) for i in range(NJ)]
MS[2][0] = var('M2l')
m = var('M0:{0}'.format(NJ))
GAM = [sympify(0) for i in range(NJ)]
J = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                            'XY{0},YY{0},YZ{0},'
                            'XZ{0},YZ{0},M{0}ll').format(i))) for i in range(NJ)]

sydi = {}
tridi = {}
for j in range(NJ):
    antTj = mat_trig_replace(transform(theta[j],r[j],alpha[j],d[j]),tridi,[j],theta,alpha,gamma)
    antAj[j] = mat_sym_replace(get_r(antTj),sydi,'A',[j])
    antPj[j] = mat_sym_replace(get_p(antTj),sydi,'L',[j])
G = Matrix([0,0,-var('G')])
w[0] = zeros(3,1)#Matrix(var('w01:4'))
wp[0] = zeros(3,1)#Matrix(var('wp01:4'))
v[0] = zeros(3,1)#Matrix(var('v01:4'))
vp[0] = -G#Matrix(var('vp01:4')) + G
for j in range(NJ):
    jAant = antAj[j].T
    qpj = Matrix([0,0,qp[j]])
    qdpj = Matrix([0,0,qdp[j]])
    if j != 0:
        wi[j] = mat_sym_replace(jAant*w[ant[j]],sydi,'WI',[j])
        w[j] = mat_sym_replace(wi[j]+(1-sigm[j])*qpj,sydi,'W',[j])
        wp[j] =  mat_sym_replace(jAant*wp[ant[j]]+(1-sigm[j])*(qdpj+hat(wi[j])*qpj),
                                sydi,'WP',[j])
        U[j] = mat_sym_replace(hat(wp[j]) + hat(w[j])*hat(w[j]),sydi,'U',[j])
        vsp[j] = mat_sym_replace(vp[ant[j]] + U[ant[j]]*antPj[j],sydi,'VSP',[j])
        vp[j] = mat_sym_replace(jAant*vsp[j] + sigm[j]*(qdpj+2*hat(wi[j])*qpj),
                                    sydi,'VP',[j])
    F[j] = mat_sym_replace(m[j]*vp[j] + U[j]*MS[j],sydi,'F',[j])
    Psi[j] = mat_sym_replace(hat(w[j])*J[j],sydi,'PSI',[j])
    N[j] = mat_sym_replace(J[j]*wp[j] + Psi[j]*wp[j] + hat(MS[j])*vp[j],
                        sydi,'No',[j])
for j in reversed(range(NJ)):
    fj[j] = mat_sym_replace(F[j] + fej[j],sydi,'E',[j])
    nj[j] = mat_sym_replace(N[j] + nej[j],sydi,'N',[j])
    i = ant[j]
    fi[j] = mat_sym_replace(antAj[j]*fj[j],sydi,'FDI',[j])
    if i != -1:
        fej[i] += fi[j]
        nej[i] += antAj[j]*nj[j] + hat(antPj[j])*fi[j]
for j in range(NJ):
    if sigm[j] != 2:
        GAM[j] = sym_replace((sigm[j]*fj[j] + (1-sigm[j])*nj[j])[2],sydi,'GAM',[j])
# + FS[j]*sign(qp[j]) +                               FV[j]*qp[j] + IA[j]*qdp[j]
for j in range(NJ):
    print  'GAM{0}'.format(j), '=', simplify(unfold(GAM[j],sydi))

