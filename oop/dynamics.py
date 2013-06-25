from sympy import sign,eye,zeros,Matrix,Integer
from symoro_oop import RX90,CartPole,trig_replace,sym_replace,sym_mat_replace,\
                         hat,get_r,get_p,sym_mats_replace
from copy import copy

chars = 'ABCDEFGHJKLMNPQRSTUVWXYZ'

def inverse_dynamic_NE(robo):
    f_ext = copy(robo.f_ext)
    n_ext = copy(robo.n_ext)
    w = init_w(robo)
    wdot,vdot = init_wv_dot(robo)
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    U = init_mat(robo)
    U.append(zeros(3,3))
    F = init_vec(robo)
    N = init_vec(robo)
    fj = init_vec(robo)
    nj = init_vec(robo)
    sydi = {}
    for j in range(robo.NL):
        compute_transform(robo,j,sydi,antRj,antPj)
    for j in range(robo.NL):
        compute_twist(robo,j,sydi,antRj,antPj,w,wdot,U,vdot)
        compute_wrench(robo,j,sydi,w,wdot,U,vdot,F,N)
    for j in reversed(range(robo.NL)):
        compute_joint_wrench(robo,j,sydi,antRj,antPj,vdot,fj,nj,F,N,f_ext,n_ext)
    for j in range(robo.NL):
        compute_torque(robo,j,sydi,fj,nj)
    return sydi

def dynamic_identification_NE(robo):
    w = init_w(robo)
    wdot,vdot = init_wv_dot(robo)
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    U = init_mat(robo)
    U.append(zeros(3,3))
    fj = init_vec(robo)
    nj = init_vec(robo)
    sydi = {}   
    robo_tmp = copy(robo)
    robo_tmp.IA = zeros(robo.NL)
    robo_tmp.FV = zeros(robo.NL)
    robo_tmp.FS = zeros(robo.NL)
    for j in range(robo.NL):
        compute_transform(robo,j,sydi,antRj,antPj)        
    for j in range(robo.NL):
        compute_twist(robo,j,sydi,antRj,antPj,w,wdot,U,vdot)
    for k in range(robo.NL):
        param_vec = get_dynam_param(robo,k)
        F = init_vec(robo)
        N = init_vec(robo)
        for i in range(10):
            if param_vec[i] == Integer(0):
                continue
            robo_tmp.num = [str(robo.num[l])+str(param_vec[i])
                            for l in range(k+1)]
            mask = zeros(10,1)
            mask[i] = 1
            put_dyn_param(robo_tmp,mask,k)
            compute_wrench(robo_tmp,k,sydi,w,wdot,U,vdot,F,N)            
            f_ext = copy(robo.f_ext)
            n_ext = copy(robo.n_ext)            
            for j in reversed(range(k+1)):
                compute_joint_wrench(robo_tmp,j,sydi,antRj,antPj,vdot,fj,nj,F,N,
                                 f_ext,n_ext)
            for j in range(k+1):
                compute_torque(robo_tmp,j,sydi,fj,nj,'DG')
        compute_joint_torque_deriv(robo.IA[k],robo.qddot[k],robo.num[k],sydi)
        compute_joint_torque_deriv(robo.FS[k],sign(robo.qdot[k]),robo.num[k],sydi)
        compute_joint_torque_deriv(robo.FV[k],robo.qdot[k],robo.num[k],sydi)
    return sydi

def direct_dynamic_NE(robo):
    wi = init_vec(robo)
    w = init_w(robo)
    jaj = init_vec(robo,6)
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    jTant = init_mat(robo,6)
    beta_star = init_vec(robo,6)
    grandJ = init_mat(robo,6)
    link_acc = init_vec(robo,6)
    H_inv = init_scalar(robo)
    juj= init_vec(robo,6) #Jj*aj/Hj
    Tau = init_scalar(robo)
    grandVp = init_vec(robo,6) 
    grandVp.append(Matrix([robo.VP0-robo.G,robo.W0]))
    sydi = {}
    for j in range(robo.NL):
        compute_transform(robo,j,sydi,antRj,antPj)
        compute_omega(robo,j,sydi,antRj,w,wi)
        compute_screw_transform(robo,j,sydi,antRj,antPj,jTant)              
        if robo.sigma[j] == 0:
            jaj[j] = Matrix([0,0,0,0,0,1])
        elif robo.sigma[j] == 1:
            jaj[j] = Matrix([0,0,1,0,0,0])
    for j in range(robo.NL):
        compute_beta(robo,j,sydi,w,beta_star)
        compute_link_acc(robo,j,sydi,antRj,antPj,link_acc,w,wi)
        grandJ[j] = inertia_spatial(robo.J[j],robo.MS[j],robo.M[j])

    for j in reversed(range(robo.NL)):
        replace_beta_J_star(robo,j,sydi,grandJ,beta_star)
        compute_Tau(robo,j,sydi,grandJ,beta_star,jaj,juj,H_inv,Tau)    
        if robo.ant[j] != -1:
            compute_beta_J_star(robo,j,sydi,grandJ,jaj,juj,Tau,
                        beta_star,jTant,link_acc)
    for j in range(robo.NL):
        compute_acceleration(robo,j,sydi,jTant,grandVp,juj,H_inv,jaj,Tau,link_acc)    
    for j in range(robo.NL):
        compute_coupled_forces(robo,j,sydi,grandVp,grandJ,beta_star)

def inertia_matrix(robo):
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    Jplus,MSplus,Mplus = init_Jplus(robo)
    AJE1 = init_vec(robo)
    f = init_vec_ext(robo)
    n = init_vec_ext(robo)
    A = zeros(robo.NL)
    sydi = {}
    for j in range(robo.NL):
        compute_transform(robo,j,sydi,antRj,antPj)
    for j in reversed(range(-1,robo.NL)):
        replace_Jplus(robo,j,sydi,Jplus,MSplus,Mplus)
        if j != -1:
            compute_Jplus(robo,j,sydi,antRj,antPj,Jplus,MSplus,Mplus,AJE1)       
    for j in range(robo.NL):      
        compute_A_diagonal(robo,j,sydi,Jplus,MSplus,Mplus,f,n,A)
        ka = j 
        while ka != -1:
            k = ka
            ka = robo.ant[ka]
            compute_A_triangle(robo,j,k,ka,sydi,antRj,antPj,f,n,A,AJE1)
    sym_mats_replace(A,sydi,'A',forced = True)
    J_base = inertia_spatial(Jplus[-1],MSplus[-1],Mplus[-1])
    sym_mats_replace(J_base,sydi,'JP',0,forced = True)            

def pseudo_force_NE(robo):
    robo_pseudo = copy(robo)
    robo_pseudo.qddot = zeros(robo_pseudo.NL)
    return inverse_dynamic_NE(robo_pseudo)

def init_U(robo):
    U = init_mat(robo)
    U.append(zeros(3,3))
    return U

def init_Jplus(robo):
    Jplus = copy(robo.J)
    Jplus.append(zeros(3))
    MSplus = copy(robo.MS)
    MSplus.append(zeros(3,1))
    Mplus = copy(robo.M)
    Mplus.append(0)
    return Jplus,MSplus,Mplus

def init_mat(robo,N = 3):
    return [zeros(N,N) for i in range(robo.NL)]
    
def init_vec(robo,N = 3):
    return [zeros(N,1) for i in range(robo.NL)]

def init_vec_ext(robo,N = 3):
    return [zeros(N,1) for i in range(robo.NL+1)]

def init_scalar(robo):
    return [0 for i in range(robo.NL)]

def init_w(robo):
    w = [zeros(3,1) for i in range(robo.NL)]
    w.append(robo.W0)
    return w
   
def init_wv_dot(robo):
    wdot = [zeros(3,1) for i in range(robo.NL)]
    wdot.append(robo.WP0)
    vdot = [zeros(3,1) for i in range(robo.NL)]
    vdot.append(robo.VP0-robo.G)
    return wdot,vdot
    
def compute_transform(robo,j,sydi,antRj,antPj):
    robo.num[j] = robo.num[j]
    antTj = robo.transform(j)
    for angle in robo.get_angles(j):
        antTj = trig_replace(antTj,sydi,angle,robo.num[j])
    antRj[j] = sym_mat_replace(get_r(antTj),sydi,'A',robo.num[j])
    antPj[j] = sym_mat_replace(get_p(antTj),sydi,'L',robo.num[j])

    
def compute_screw_transform(robo,j,sydi,antRj,antPj,jTant):
    jRant = antRj[j].T 
    ET = sym_mat_replace(-jRant*hat(antPj[j]),sydi,'JPR',robo.num[j])
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3).row_join(jRant)]))
                        
def compute_twist(robo,j,sydi,antRj,antPj,w,wdot,U,vdot):
    jRant = antRj[j].T
    qdj = Matrix([0,0,robo.qdot[j]])
    qddj = Matrix([0,0,robo.qddot[j]])
    wi = sym_mat_replace(jRant*w[robo.ant[j]],sydi,'WI',robo.num[j])
    w[j] = sym_mat_replace(wi+(1-robo.sigma[j])*qdj,sydi,'W',robo.num[j])
    wdot[j] = jRant*wdot[robo.ant[j]]+(1-robo.sigma[j])*(qddj+hat(wi)*qdj)
    sym_mat_replace(wdot[j],sydi,'WP',robo.num[j])                   
    U[j] = hat(w[j])*hat(w[j])+hat(wdot[j])
    sym_mat_replace(U[j],sydi,'U',robo.num[j])    
    vsp = vdot[robo.ant[j]]+U[robo.ant[j]]*antPj[j]
    sym_mat_replace(vsp,sydi,'VSP',robo.num[j])
    vdot[j] = jRant*vsp+robo.sigma[j]*(qddj+2*hat(wi)*qdj)
    sym_mat_replace(vdot[j],sydi,'VP',robo.num[j])  
                                
def compute_omega(robo,j,sydi,antRj,w,wi):
    jRant = antRj[j].T
    qdj = Matrix([0,0,robo.qdot[j]])
    wi[j] = sym_mat_replace(jRant*w[robo.ant[j]],sydi,'WI',robo.num[j])
    w[j] = sym_mat_replace(wi[j]+(1-robo.sigma[j])*qdj,sydi,'W',robo.num[j])
                             
def compute_wrench(robo,j,sydi,w,wdot,U,vdot,F,N):
    F[j] = robo.M[j]*vdot[j]+U[j]*robo.MS[j]
    sym_mat_replace(F[j],sydi,'F',robo.num[j])
    Psi = sym_mat_replace(robo.J[j]*w[j],sydi,'PSI',robo.num[j])
    N[j] = robo.J[j]*wdot[j]+hat(w[j])*Psi
    sym_mat_replace(N[j],sydi,'No',robo.num[j])
                        

def compute_joint_wrench(robo,j,sydi,antRj,antPj,vdot,fj,nj,F,N,f_ext,n_ext):
    fj[j] = sym_mat_replace(F[j]+f_ext[j],sydi,'E',robo.num[j])
    nj[j] = N[j]+n_ext[j]+hat(robo.MS[j])*vdot[j]
    sym_mat_replace(nj[j],sydi,'N',robo.num[j])
    f_ant = sym_mat_replace(antRj[j]*fj[j],sydi,'FDI',robo.num[j])
    if robo.ant[j] != -1:
        f_ext[robo.ant[j]] += f_ant
        n_ext[robo.ant[j]] += antRj[j]*nj[j]+hat(antPj[j])*f_ant

def compute_torque(robo,j,sydi,fj,nj,name = 'GAM'):
    if robo.sigma[j] != 2:
        tau = (robo.sigma[j]*fj[j]+(1-robo.sigma[j])*nj[j])
        tau_total = tau[2]+robo.fric_s(j)+robo.fric_v(j)+robo.tau_ia(j)
        tau_total = sym_replace(tau_total,sydi,name,robo.num[j],forced = True)
  
def put_dyn_param(robo,P,j):
    robo.J[j] = Matrix([[P[0],P[1],P[2]],
                [P[1],P[3],P[4]],
                [P[2],P[4],P[5]]])
    robo.MS[j] = Matrix(3,1,P[6:9])
    robo.M[j] = P[9]

def get_dynam_param(robo,j):
    return [robo.J[j][0],robo.J[j][1],robo.J[j][2],robo.J[j][4],
                    robo.J[j][5],robo.J[j][8],robo.MS[j][0],robo.MS[j][1],
                    robo.MS[j][2],robo.M[j]]

def inertia_spatial(J,MS,M):
    return Matrix([(M*eye(3)).row_join(hat(MS).T),hat(MS).row_join(J)])

def compute_joint_torque_deriv(param,arg,index,sydi):
    if param != Integer(0) and arg != Integer(0):
        index = str(index)+str(param)
        arg = sym_replace(arg,sydi,'DG',index,forced=True)

def compute_beta(robo,j,sydi,w,beta_star):
    E1 = sym_mat_replace(robo.J[j]*w[j],sydi,'JW',robo.num[j])
    E2 = sym_mat_replace(hat(w[j])*E1,sydi,'KW',robo.num[j])
    E3 = hat(w[j])*robo.MS[j]
    E4 = sym_mat_replace(hat(w[j])*E3,sydi,'SW',robo.num[j])
    E5 = -robo.n_ext[j]-E2
    E6 = -robo.f_ext[j]-E4
    beta_star[j] = Matrix([E6,E5])

def compute_link_acc(robo,j,sydi,antRj,antPj,link_acc,w,wi):
    E1 = sym_mat_replace(hat(wi[j])*Matrix([0,0,robo.qdot[j]]),
                             sydi,'WQ',robo.num[j])
    E2 = (1-robo.sigma[j])*E1
    E3 = 2*robo.sigma[j]*E1
    E4 = hat(w[robo.ant[j]])*antPj[j]
    E5 = hat(w[robo.ant[j]])*E4
    E6 = antRj[j].T*E5
    E7 = sym_mat_replace(E6+E3,sydi,'LW',robo.num[j])
    link_acc[j] = Matrix([E7,E2])

def replace_beta_J_star(robo,j,sydi,grandJ,beta_star):
    grandJ[j] = sym_mats_replace(grandJ[j],sydi,'MJE',robo.num[j])
    beta_star[j] = sym_mat_replace(beta_star[j],sydi,'VBE',robo.num[j])
  
def compute_Tau(robo,j,sydi,grandJ,beta_star,jaj,juj,H_inv,Tau):
    Jstar_jaj = grandJ[j]*jaj[j]
    if robo.sigma[j] == 2:
        Tau[j] = 0
    else:
        H_inv[j] = 1/(jaj[j].dot(Jstar_jaj)+robo.IA[j])
        H_inv[j] = sym_replace(H_inv[j],sydi,'JD',robo.num[j])
        juj[j] = Jstar_jaj*H_inv[j]
        sym_mat_replace(juj[j],sydi,'JU',robo.num[j])
        joint_friction = robo.fric_s(j)+robo.fric_v(j)
        Tau[j] = jaj[j].dot(beta_star[j])+robo.GAM[j]-joint_friction
        Tau[j] = sym_replace(Tau[j],sydi,'GW',robo.num[j])

def compute_beta_J_star(robo,j,sydi,grandJ,jaj,juj,Tau,
                        beta_star,jTant,link_acc):
    Jstar_jaj = grandJ[j]*jaj[j]
    grandK = sym_mat_replace(grandJ[j]-juj[j]*Jstar_jaj.T,sydi,'GK',robo.num[j])
    E1 = sym_mat_replace(grandK*link_acc[j],sydi,'NG',robo.num[j])
    E3 = sym_mat_replace(E1 + Tau[j]*juj[j],sydi,'VS',robo.num[j])
    alpha = sym_mat_replace(E3 - beta_star[j],sydi,'AP',robo.num[j])
    E4 = sym_mat_replace(jTant[j].T*grandK,sydi,'GX',robo.num[j])
    E5 = sym_mats_replace(E4*jTant[j],sydi,'TKT',robo.num[j])
    grandJ[robo.ant[j]] += E5
    beta_star[robo.ant[j]] -= jTant[j].T*alpha
   
def compute_acceleration(robo,j,sydi,jTant,grandVp,juj,H_inv,jaj,Tau,link_acc):
    grandR = sym_mat_replace(jTant[j]*grandVp[robo.ant[j]]+link_acc[j],
                              sydi,'VR',robo.num[j])       
    E1 = sym_replace(juj[j].dot(grandR),sydi,'GU',robo.num[j])
    if robo.sigma[j] == 2: 
        qddot = 0
    else:
        qddot = H_inv[j]*Tau[j]- E1       
    qddot = sym_replace(qddot,sydi,"QDP",robo.num[j],forced = True)
    grandVp[j] = (grandR + qddot*jaj[j]) 
    grandVp[j][3:,0] = sym_mat_replace(grandVp[j][3:,0],sydi,'WP',robo.num[j])
    grandVp[j][:3,0] = sym_mat_replace(grandVp[j][:3,0],sydi,'VP',robo.num[j])
    
def compute_coupled_forces(robo,j,sydi,grandVp,grandJ,beta_star):
    E3 = sym_mat_replace(grandJ[j]*grandVp[j],sydi,'DY',robo.num[j])
    couplforce = E3 - beta_star[j];
    sym_mat_replace(couplforce[3:,0],sydi,'N',robo.num[j])
    sym_mat_replace(couplforce[:3,0],sydi,'E',robo.num[j])
    
def replace_Jplus(robo,j,sydi,Jplus,MSplus,Mplus):
    sym_mat_replace(Jplus[j],sydi,'JP',robo.num[j])
    sym_mat_replace(MSplus[j],sydi,'MSP',robo.num[j])
    Mplus[j] = sym_replace(Mplus[j],sydi,'MP',robo.num[j])

def compute_Jplus(robo,j,sydi,antRj,antPj,Jplus,MSplus,Mplus,AJE1):
    hat_antPj = hat(antPj[j])
    antMSj = sym_mat_replace(antRj[j]*MSplus[j],sydi,'AS',robo.num[j])
    E1 = sym_mat_replace(antRj[j]*Jplus[j],sydi,'AJ',robo.num[j])
    AJE1[j] = E1[:,2]
    E2 = sym_mat_replace(E1*antRj[j].T,sydi,'AJA',robo.num[j])
    E3 = sym_mat_replace(hat_antPj*hat(antMSj),sydi,'PAS',robo.num[j])    
    Jplus[robo.ant[j]] += E2-(E3+E3.T)+hat_antPj*hat_antPj.T*Mplus[j]
    MSplus[robo.ant[j]] += antMSj+antPj[j]*Mplus[j]
    Mplus[robo.ant[j]] += Mplus[j]

def compute_A_diagonal(robo,j,sydi,Jplus,MSplus,Mplus,f,n,A):
    if robo.sigma[j]==0:
        f[j]=Matrix([-MSplus[j][1],MSplus[j][0],0])                                 
        n[j]=Jplus[j][:,2]                       
        A[j,j]=Jplus[j][2,2]+robo.IA[j]        
    elif robo.sigma[j] == 1:
        f[j]=Matrix([0,0,Mplus[j]])
        n[j]=Matrix([MSplus[j][1],-MSplus[j][0],0])
        A[j,j]=Mplus[j]+robo.IA[j]  
    sym_mat_replace(f[j],sydi,'E'+chars[j],robo.num[j])
    sym_mat_replace(n[j],sydi,'N'+chars[j],robo.num[j])
     
def compute_A_triangle(robo,j,k,ka,sydi,antRj,antPj,f,n,A,AJE1):
    f[ka] = antRj[k]*f[k]
    if k == j and robo.sigma[j] ==0:
        n[ka] = AJE1[k]+hat(antPj[k])*f[k]
    else:
        n[ka] = antRj[k]*n[k]+hat(antPj[k])*f[k]
    if ka == -1:
        sym_mat_replace(f[ka],sydi,'AV0')
        sym_mat_replace(n[ka],sydi,'AW0')
    else:
        sym_mat_replace(f[ka],sydi,'E'+chars[j],robo.num[ka])
        sym_mat_replace(n[ka],sydi,'N'+chars[j],robo.num[ka]) 
        if robo.sigma[ka] == 0:
            A[j,ka] = n[ka][2]
        elif robo.sigma[ka] == 1:
            A[j,ka] = f[ka][2]
                    
print 'Inverse dynamic model using Newton-Euler Algorith'
inverse_dynamic_NE(CartPole())

print 'Pseudo forces using Newton-Euler Algorith'
pseudo_force_NE(CartPole())

print 'Direct dynamic model using Newton-Euler Algorith'
direct_dynamic_NE(CartPole())

print 'Dynamic identification model using Newton-Euler Algorith'
dynamic_identification_NE(CartPole())

print 'Inertia Matrix using composite links'
inertia_matrix(CartPole())

