# -*- coding: utf-8 -*-
"""
Created on Wed May 29 20:26:47 2013

@author: Bogdqn
"""
from sympy import sign, zeros
from symoro import *
from acrobat import *

wi = [zeros(3,1) for i in range(NL+1)]
w = [zeros(3,1) for i in range(NL+1)]
phij = [zeros(6,1) for i in range(NL+1)]
lamdaj = [0 for i in range(NL+1)]
effejppj = [zeros(3,1) for i in range(NL+1)]
antAj = [zeros(3,3) for i in range(NL+1)]
antPj = [zeros(3,1) for i in range(NL+1)]
jTant = [zeros(6,6) for i in range(NL+1)]

sydi = {}
tridi = {}
for j in range(NL):
    index = num[j]
    antTj = mat_trig_replace(transform(theta[j],r[j],alpha[j],d[j]),tridi,[j],theta,alpha,gamma)
    antAj[j] = mat_sym_replace(get_r(antTj),sydi,'A',index)
    antPj[j] = mat_sym_replace(get_p(antTj),sydi,'L',index)
    jAant = antAj[j].T
    
    if ant[j] != -1:
        want = w[ant[j]]
    else:
        want = W0
        
    wi[j] = mat_sym_replace(jAant*want,sydi,'WI',index)
    w[j] = mat_sym_replace(wi[j] + (1-sigm[j])*Matrix([0,0,qp[j]]),sydi,'W',index)
    quantT = mat_sym_replace(-jAant*hat(antPj[j]),sydi,'JPR',index)
    jTant[j] = (Matrix([jAant.row_join(zeros(3)),
                    quantT.row_join(jAant)]))
    if sigm[j] == 1:
        phij[j] = Matrix([0,0,0,0,0,1])
        lamdaj[j] = GAM[j]
        effejppj[j] = Matrix([0,0,0,0,0,GAM[j]])
    elif sigm[j] == 0:
        phij[j] = Matrix([0,0,1,0,0,0])
        lamdaj[j] = GAM[j]
        effejppj[j] = Matrix([0,0,GAM[j],0,0,0])
    elif sigm[j] == 2:
        phij[j] = Matrix([0,0,0,0,0,0])
        lamdaj[j] = GAM[j]
        effejppj[j] = Matrix([0,0,0,0,0,0])

vecbetetoil = [zeros(6,1) for i in range(NL+1)]
matMJetoil = [zeros(6,6) for i in range(NL+1)]
gamaj = [zeros(6,1) for i in range(NL+1)]

for j in range(NL):
    index = num[j]
    member1 = mat_sym_replace(J[j]*w[j],sydi,'JW',index)
    member2 = mat_sym_replace(hat(w[j])*member1,sydi,'KW',index)
    member3 = hat(w[j])*MS[j]
    member4 = mat_sym_replace(hat(w[j])*member3,sydi,'SW',index)
    member5 = -nej[j]-member2
    member6 = -fej[j]-member4
    vecbetetoil[j] = Matrix([member5,member6])
    if ant[j] != -1:
        want = w[ant[j]]
    else:
        want = W0
    quant1 = mat_sym_replace(hat(wi[j])*Matrix([0,0,qp[j]]),sydi,'WQ',index)
    quant2 = (1-sigm[j]) * quant1
    quant3 = 2*sigm[j] * quant1
    quant4 = hat(want)*antPj[j]
    quant5 = hat(want)*quant4
    quant6 = antAj[j].T*quant5
    quant7 = mat_sym_replace(quant6+quant3,sydi,'LW',index)
    gamaj[j] = Matrix([quant2,quant7])

    #TODO ask about this matrix, the meaning
    matMJetoil[j] = Matrix([J[j].row_join(hat(MS[j])),hat(MS[j]).T.row_join(M[j]*eye(3))])

jdej = [0 for i in range(NL+1)] #1/Hj
juj= [zeros(6,1) for i in range(NL+1)] #Jj*aj/Hj
jvej = [zeros(6,1) for i in range(NL+1)] #Jj*aj
grandKj = [zeros(6) for i in range(NL+1)] 
grandwj = [0 for i in range(NL+1)]
alphaj = [zeros(6,1) for i in range(NL+1)]

#phi is jaj
for j in reversed(range(NL)):
    index = num[j]
    jvej[j] = matMJetoil[j]*phij[j] 
    matMJetoil[j] = matsymm_sym_replace(matMJetoil[j],sydi,'MJE',index)
    vecbetetoil[j] = mat_sym_replace(vecbetetoil[j],sydi,'VBE',index)
    if sigm[j] == 2:
        grandwj[j] = 0
    else:
        jdej[j] = sym_replace(1/(phij[j].dot(jvej[j])+IA[j]),sydi,'JD',index)
        juj[j] = mat_sym_replace(jdej[j]*jvej[j],sydi,'JU',index)
        grandwj[j] = sym_replace(phij[j].dot(vecbetetoil[j])+lamdaj[j]-FS[j]*sign(qp[j])-FV[j]*qp[j],
                   sydi,'GW',index)
    


    if ant[j] != -1:
        #    	(***************************************)
        #     (*  calculculation of"projected matrix"*)
        #     (***************************************)
        #??? why is this inside the condition section?
        grandKj[j] = mat_sym_replace(matMJetoil[j]-juj[j]*jvej[j].T,sydi,'GK',index)
        quant1 = mat_sym_replace(grandKj[j]*gamaj[j],sydi,'NG',index) 
        quant3 = mat_sym_replace(quant1 + grandwj[j]*juj[j],sydi,'VS',index)
        alphaj[j] = mat_sym_replace(quant3 - vecbetetoil[j],sydi,'AP',index)
        
        #     (***************************************)
        #      (*  calculation of"effective matrix"*)
        #     (***************************************)
        quant4 = mat_sym_replace(jTant[j].T*grandKj[j],sydi,'GX',index)
        quant5 = matsymm_sym_replace(quant4*jTant[j],sydi,'TKT',index)
        #??? maybe separated routine is needed for replacement to evoid recursion
        matMJetoil[ant[j]] = matMJetoil[ant[j]]+quant5
        vecbetetoil[ant[j]] = vecbetetoil[ant[j]]-jTant[j].T*alphaj[j]
    
qdp = [0 for i in range(NL+1)]
grandVp = [zeros(6,1) for i in range(NL+1)]
grandRj = [zeros(6,1) for i in range(NL+1)]

for j in range(NL):
    index = num[j]
    if ant[j] != -1:
        grandVpant = grandVp[ant[j]]
    else:
        grandVpant = Matrix([WP0,VP0-G])
    grandRj[j] = mat_sym_replace(jTant[j]*grandVpant+gamaj[j],sydi,'VR',index)
       
    member1 = sym_replace(juj[j].dot(grandRj[j]),sydi,'GU',index)

    if sigm[j] == 2: 
        qdp[j] = 0
    else:
        qdp[j] = jdej[j]*grandwj[j]- member1
       
    qdp[j] = sym_replace(qdp[j],sydi,"QDP",index)
    grandVp[j] = (grandRj[j] + qdp[j]*phij[j]) 
    grandVp[j][:3,0] = mat_sym_replace(grandVp[j][:3,0],sydi,'WP',index)
    grandVp[j][3:,0] = mat_sym_replace(grandVp[j][3:,0],sydi,'VP',index)
 
for j in range(NL):
    index = num[j]
    member3 = mat_sym_replace(matMJetoil[j]*grandVp[j],sydi,'DY',index)
    couplforce = member3 - vecbetetoil[j];
    mat_sym_replace(couplforce[:3,0],sydi,'N',index)
    mat_sym_replace(couplforce[3:,0],sydi,'E',index)
            
       