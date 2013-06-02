# -*- coding: utf-8 -*-
"""
Created on Wed May 29 20:26:47 2013

@author: Bogdqn
"""
from sympy import pi, sign
from symoro import *

#needed variables
D3 = var('D3')
RL4 = var('RL4')
#table of geometric parameters RX90
ant = (-1,0,1,2,3,4)
sigm = (0,0,0,0,0,0)
alpha = (0,pi/2,0,-pi/2,pi/2,-pi/2)
d = (0,0,D3,0,0,0)
theta = var('th1:7')
r = (0,0,0,RL4,0,0)
b = (0,0,0,0,0,0)
gamma = (0,0,0,0,0,0)
q = var('q1:7')
qp = var('qp1:7')
omegajant = [zeros(3,1) for i in range(6)]
omegaj = [zeros(3,1) for i in range(6)]
phij = [zeros(6,1) for i in range(6)]
lamdaj = [0 for i in range(6)]
effappj = [zeros(3,1) for i in range(6)]
antAj = [zeros(3,3) for i in range(6)]
antPj = [zeros(3,1) for i in range(6)]
jKant = [zeros(6,6) for i in range(6)]
MS = [Matrix(var('MX{0},MY{0},MZ{0}'.format(i+1))) for i in range(6)]
GAM = var('GAM1:7')   #Torques
na = [Matrix(var('CX{0},CY{0},CZ{0}'.format(i+1))) for i in range(6)]
fa = [Matrix(var('FX{0},FY{0},FZ{0}'.format(i+1))) for i in range(6)]
matJj = [Matrix(3,3,var(('XX{0},XY{0},XZ{0},'
                            'XY{0},YY{0},YZ{0},'
                            'XZ{0},YZ{0},ZZ{0}').format(i+1))) for i in range(6)]
sydi = {}
tridi = {}
for j in range(6):
    antTj = mat_trig_replace(transform(theta[j],r[j],alpha[j],d[j]),tridi,[j],theta)
    antAj[j] = mat_sym_replace(get_r(antTj),sydi,'A',[j])
    antPj[j] = mat_sym_replace(get_p(antTj),sydi,'L',[j])
    jAant = antAj[j].T
    antPjpre = hat(antPj[j])
    if j != 0:
        omegajant[j] = mat_sym_replace(jAant*omegaj[ant[j]],sydi,'WI',[j])
    omegaj[j] = mat_sym_replace(omegajant[j] + (1-sigm[j])*Matrix([0,0,qp[j]]),sydi,'W',[j])
    jKant[j] = (Matrix([jAant.row_join(-jAant*antPjpre),
                    zeros(3).row_join(jAant)]))
    if sigm[j] == 1:
        phij[j] = Matrix([0,0,0,0,0,1])
        lamdaj[j] = GAM[j]
        effappj[j] = Matrix([0,0,0,0,0,GAM[j]])
    elif sigm[j] == 0:
        phij[j] = Matrix([0,0,1,0,0,0])
        lamdaj[j] = GAM[j]
        effappj[j] = Matrix([0,0,GAM[j],0,0,0])
    elif sigm[j] == 2:
        phij[j] = Matrix([0,0,0,0,0,0])
        lamdaj[j] = GAM[j]
        effappj[j] = Matrix([0,0,0,0,0,0])

vecbetetoil = [zeros(6,1) for i in range(6)]
matMJetoil = [zeros(6,6) for i in range(6)]
gamaj = [zeros(6,1) for i in range(6)]
M = var('M1:7')#???
for j in range(6):
    member1 = mat_sym_replace(matJj[j]*omegaj[j],sydi,'JW',[j])
    member2 = mat_sym_replace(hat(omegaj[j])*member1,sydi,'KW',[j])
    member3 = hat(omegaj[j])*MS[j]
    member4 = mat_sym_replace(hat(omegaj[j])*member3,sydi,'SW',[j])
    member5 = -na[j]-member2
    member6 = -fa[j]-member4
    vecbetetoil[j] = Matrix([member5,member6])

    quant1 = mat_sym_replace(hat(omegajant[j])*Matrix([0,0,qp[j]]),sydi,'WQ',[j])
    quant2 = (1-sigm[j]) * quant1
    quant3 = 2*sigm[j] * quant1
    quant4 = hat(omegaj[ant[j]])*antPj[j]
    quant5 = hat(omegaj[ant[j]])*quant4
    quant6 = antAj[j].T*quant5
    quant7 = mat_sym_replace(quant6+quant3,sydi,'LW',[j])
    gamaj[j] = Matrix([quant2,quant7])

    MSjpre=hat(MS[j])
    quant8=MSjpre.T
    #TODO ask about this matrix, the meaning
    matMJetoil[j] = Matrix([matJj[j].row_join(quant8),MSjpre.row_join(M[j]*eye(3))])

jdej = [0 for i in range(6)]
juj = [zeros(6,1) for i in range(6)]
jvej = [zeros(6,1) for i in range(6)]
grandwj = [0 for i in range(6)]
IA = var('IA1:7')
FS = var('FS1:7')
FV = var('FV1:7')
for j in reversed(range(6)):
    jvej[j] = matMJetoil[j]*phij[j]
    if sigm[j]==2:
        grandwj[j]=0
    else:
        jdej[j] = sym_replace(1/(phij[j].dot(jvej[j])+IA[j]),sydi,'JD',[j])#????
        juj[j] = mat_sym_replace(jdej[j] * jvej[j],sydi,'JU',[j])
        grandwj[j] = sym_replace(phij[j].dot(vecbetetoil[j])+lamdaj[j]-FS[j]*sign(qp[j])-FV[j]*qp[j],
                   sydi,'GW',[j])



#    if ant[j]!=0:
#    	(***************************************)
#     (*  calculculation of"projected matrix"*)
#     (***************************************)
#       expr1={{juj[[e1,1]]},{juj[[e1,2]]},{juj[[e1,3]]},
#      {juj[[e1,4]]},{juj[[e1,5]]},{juj[[e1,6]]}};
#
#       expr= expr1.{jvej[j]};
#
#       grandNj[j]=matMJetoil[j]-expr;
#     grandNj[j]=ecritmatrice[grandNj[j],"GK",3,e1,6,6];
#       quant1= grandNj[j].gamaj[j];
#
#       quant1=ecritvecteur[quant1,"NG",3,e1,6];
#       quant2= grandwj[j] juj[j];
#       quant3= quant1+quant2;
#
#       quant3= ecritvecteur[quant3,"VS",3,e1,6];
#       alphaj[j]=quant3 - vecbetetoil[j];
#
#       alphaj[j]=ecritvecteur[ alphaj[j],"AP",3,e1,6];
#
#
#     (***************************************)
#      (*  calculation of"effective matrix"*)
#     (***************************************)
#
#     quant4=Transpose[jxant[j]].grandNj[j];
#     quant4=ecritmatrice[quant4,"GX",3,e1,6,6];
#     quant5= quant4.jxant[j];
#     quant5=ecrituresymetrie[quant5,"TKT",3,e1,6];
#
#    matMJetoil[[Ant[j]]]=matMJetoil[[Ant[j]]]+quant5;
#    matMJetoil[[Ant[j]]]=
#    ecrituresymetrie[matMJetoil[[Ant[j]]],"MJE",3,Ant[j],6];
#
#
#
#
#
#   vecbetetoil[[Ant[j]]]=vecbetetoil[[Ant[j]]]-
#   Transpose[jxant[j]].alphaj[j];
#  vecbetetoil[[Ant[j]]]=ecritvecteur[vecbetetoil[[Ant[j]]],"VBE",3,Ant[j],6]
