"""
This module of SYMORO package provides symbolic 
modeling of robots' dynamics.

The core symbolic library is sympy.
Needed modules : symoro_oop.py

ECN - ARIA1 2013
"""
from sympy import sign,eye,zeros,Matrix,Integer
from symoro import Symoro,RX90,CartPole
from copy import copy
from geometric import transform,dgm_serial

chars = 'ABCDEFGHJKLMNPQRSTUVWXYZ'

def Newton_Euler(robo,symo):
    """Internal function. Computes Inverse Dynamic Model using
    Newton-Euler formulation
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    symo : Symoro
        Instance of symbolic manager
    """
    #init external forces
    Fex = copy(robo.Fex) 
    Nex = copy(robo.Nex)
    #init velocities and accelerations
    w = init_w(robo)
    wdot,vdot = init_wv_dot(robo)
    #init transformation 
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    #init auxilary matrix
    U = init_U(robo) 
    #init forces vectors
    F = init_vec(robo)
    N = init_vec(robo)
    Fjnt = init_vec(robo)
    Njnt = init_vec(robo)
    for j in range(robo.NL):
        compute_transform(robo,symo,j,antRj,antPj)
    for j in range(robo.NL):
        compute_twist(robo,symo,j,antRj,antPj,w,wdot,U,vdot)
        compute_wrench(robo,symo,j,w,wdot,U,vdot,F,N)
    for j in reversed(range(robo.NL)):
        compute_joint_wrench(robo,symo,j,antRj,antPj,vdot,Fjnt,Njnt,F,N,Fex,Nex)
    for j in range(robo.NL):
        compute_torque(robo,symo,j,Fjnt,Njnt)

def dynamic_identification_NE(robo):
    """Computes Dynamic Identification Model using
    Newton-Euler formulation
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    #init velocities and accelerations
    w = init_w(robo)
    wdot,vdot = init_wv_dot(robo)
    #init transformation 
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    #init auxilary matrix (3x3)
    U = init_U(robo)
    #init forces vectors
    Fjnt = init_vec(robo)
    Njnt = init_vec(robo)
    #init file output, writing the robot description
    symo = Symoro()   
    symo.file_open(robo,'dim')
    title = "Dynamic identification model using Newton - Euler Algorith"
    symo.write_dynam_param(robo,title)
    #virtual robot with only one non-zero parameter at once
    robo_tmp = copy(robo)
    robo_tmp.IA = zeros(robo.NL,1)
    robo_tmp.FV = zeros(robo.NL,1)
    robo_tmp.FS = zeros(robo.NL,1)
    for j in range(robo.NL):
        compute_transform(robo,symo,j,antRj,antPj)        
    for j in range(robo.NL):
        compute_twist(robo,symo,j,antRj,antPj,w,wdot,U,vdot)
    for k in range(robo.NL):
        param_vec = robo.get_dynam_param(k)
        F = init_vec(robo)
        N = init_vec(robo)
        for i in range(10):
            if param_vec[i] == Integer(0):
                continue
            #change link names according to current non-zero parameter
            robo_tmp.num = [str(robo.num[l]) + str(param_vec[i])
                            for l in range(k + 1)]
            #set the parameter to 1
            mask = zeros(10,1)
            mask[i] = 1
            robo_tmp.put_dynam_param(mask,k)
            #compute the total forcec of the link k
            compute_wrench(robo_tmp,symo,k,w,wdot,U,vdot,F,N)
            #init external forces
            Fex = copy(robo.Fex)
            Nex = copy(robo.Nex)            
            for j in reversed(range(k + 1)):
                compute_joint_wrench(robo_tmp,symo,j,antRj,antPj,
                                     vdot,Fjnt,Njnt,F,N,Fex,Nex)
            for j in range(k + 1):
                compute_torque(robo_tmp,symo,j,Fjnt,Njnt,'DG')
        #reset all the parameters to zero
        robo_tmp.put_dynam_param(zeros(10,1),k)
        #compute model for the joint parameters
        compute_joint_torque_deriv(symo,robo.IA[k],robo.qddot[k],robo.num[k])
        compute_joint_torque_deriv(symo,robo.FS[k],sign(robo.qdot[k]),robo.num[k])
        compute_joint_torque_deriv(symo,robo.FV[k],robo.qdot[k],robo.num[k])
    #closing the output file
    symo.file_out.close()
    return symo.sydi

def direct_dynamic_NE(robo):
    """Computes Direct Dynamic Model using
    Newton-Euler formulation
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    wi = init_vec(robo) #antecedent angular velocity, projected into jth frame 
    w = init_w(robo)
    jaj = init_vec(robo,6)
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    jTant = init_mat(robo,6) #Twist transform list of Matrices 6x6  
    beta_star = init_vec(robo,6)
    grandJ = init_mat(robo,6)
    link_acc = init_vec(robo,6)
    H_inv = init_scalar(robo)
    juj= init_vec(robo,6) #Jj*aj / Hj
    Tau = init_scalar(robo)
    grandVp = init_vec(robo,6) 
    grandVp.append(Matrix([robo.vdot0 - robo.G,robo.w0]))
    symo = Symoro()
    symo.file_open(robo,'ddm')
    symo.write_dynam_param(robo,
                           'Direct dynamic model using Newton - Euler Algorith')
    for j in range(robo.NL):
        compute_transform(robo,symo,j,antRj,antPj)
        compute_omega(robo,symo,j,antRj,w,wi)
        compute_screw_transform(robo,symo,j,antRj,antPj,jTant)
        if robo.sigma[j] == 0:
            jaj[j] = Matrix([0,0,0,0,0,1])
        elif robo.sigma[j] == 1:
            jaj[j] = Matrix([0,0,1,0,0,0])
    for j in range(robo.NL):
        compute_beta(robo,symo,j,w,beta_star)
        compute_link_acc(robo,symo,j,antRj,antPj,link_acc,w,wi)
        grandJ[j] = inertia_spatial(robo.J[j],robo.MS[j],robo.M[j])
    for j in reversed(range(robo.NL)):
        replace_beta_J_star(robo,symo,j,grandJ,beta_star)
        compute_Tau(robo,symo,j,grandJ,beta_star,jaj,juj,H_inv,Tau)    
        if robo.ant[j] != - 1:
            compute_beta_J_star(robo,symo,j,grandJ,jaj,juj,Tau,
                        beta_star,jTant,link_acc)
    for j in range(robo.NL):
        compute_acceleration(robo,symo,j,jTant,grandVp,juj,H_inv,jaj,Tau,link_acc) 
    for j in range(robo.NL):
        compute_coupled_forces(robo,symo,j,grandVp,grandJ,beta_star)
    symo.file_out.close()
    return symo.sydi

def inertia_matrix(robo):
    """Computes Inertia Matrix using composed link
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    Jplus,MSplus,Mplus = init_Jplus(robo)
    AJE1 = init_vec(robo)
    f = init_vec(robo,ext=1)
    n = init_vec(robo,ext=1)
    A = zeros(robo.NL,robo.NL)
    symo = Symoro()
    symo.file_open(robo,'inm')
    symo.write_dynam_param(robo,'Inertia Matrix using composite links')
    for j in range(robo.NL):
        compute_transform(robo,symo,j,antRj,antPj)
    for j in reversed(range( - 1,robo.NL)):
        replace_Jplus(robo,symo,j,Jplus,MSplus,Mplus)
        if j != - 1:
            compute_Jplus(robo,symo,j,antRj,antPj,Jplus,MSplus,Mplus,AJE1)
    for j in range(robo.NL):      
        compute_A_diagonal(robo,symo,j,Jplus,MSplus,Mplus,f,n,A)
        ka = j 
        while ka != - 1:
            k = ka
            ka = robo.ant[ka]
            compute_A_triangle(robo,symo,j,k,ka,antRj,antPj,f,n,A,AJE1)
    symo.mat_replace(A,'A',forced = True,symmet = True)
    J_base = inertia_spatial(Jplus[ - 1],MSplus[ - 1],Mplus[ - 1])
    symo.mat_replace(J_base,'JP',0,forced = True,symmet = True)      
    symo.file_out.close()
    return symo.sydi      

def inverse_dynamic_NE(robo):
    """Computes Inverse Dynamic Model using
    Newton-Euler formulation
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    symo = Symoro()
    symo.file_open(robo,'idm')
    symo.write_dynam_param(robo,
                           'Inverse dynamic model using Newton - Euler Algorith')
    Newton_Euler(robo,symo)
    symo.file_out.close()
    return symo.sydi

def pseudo_force_NE(robo):
    """Computes Coriolis, Centrifugal, Gravity, Friction and external
    torques using Newton-Euler formulation
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    robo_pseudo = copy(robo)
    robo_pseudo.qddot = zeros(robo_pseudo.NL,1)
    symo = Symoro()
    symo.file_open(robo,'ccg')
    symo.write_dynam_param(robo,
                           'Pseudo forces using Newton - Euler Algorith')
    Newton_Euler(robo_pseudo,symo)
    symo.file_out.close()
    return symo.sydi

def get_r(T):
    """Extracts rotation 3x3 matrix from 4x4 transformation matrix
    
    Parameters
    ==========
    T : Matrix 4x4
        Transformation matrix
        
    Returns
    =======
    get_r : Matrix 3x3
    """
    return T[:3,:3]

def get_p(T):
    """Extracts translation vector from 4x4 transformation matrix
    
    Parameters
    ==========
    T : Matrix 4x4
        Transformation matrix
        
    Returns
    =======
    get_p : Matrix 3x1
    """
    return T[:3,3]

def hat(v):
    """Generates vectorial preproduct matrix
    
    Parameters
    ==========
    v : Matrix 3x1
        vector
        
    Returns
    =======
    hat : Matrix 3x3
    """
    return Matrix([[0, - v[2],v[1]],
                   [v[2],0, - v[0]],
                   [ - v[1],v[0],0]])

def product_combinations(v):
    """Generates 6-vector of different v elements' product combinations
    
    Parameters
    ==========
    v : Matrix 3x1
        vector
        
    Returns
    =======
    product_combinations : Matrix 6x1
    """
    return Matrix([v[0]*v[0],v[0]*v[1],v[0]*v[2],
                 v[1]*v[1],v[1]*v[2],v[2]*v[2]]) 
                 
def init_U(robo):
    """Generates a list of auxiliary U matrices"""
    U = init_mat(robo)
    #the value for the -1th base frame
    U.append(hat(robo.w0)**2 + hat(robo.wdot0))  
    return U

def init_Jplus(robo):
    """Copies the inertia parameters. Used for composed link inertia computation

    Returns
    =======
    Jplus: list of Matrices 3x3
    MSplus: list of Matrices 3x1
    Mplus : list of var
    """
    Jplus = copy(robo.J)
    Jplus.append(zeros(3,3))
    MSplus = copy(robo.MS)
    MSplus.append(zeros(3,1))
    Mplus = copy(robo.M)
    Mplus.append(0)
    return Jplus,MSplus,Mplus

def init_mat(robo,N = 3):
    """Generates a list of Matrices.Size of the list is number of links.
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    N : int, optional
        size of the matries, default is 3
        
    Returns
    =======
    list of Matrices NxN
    """
    return [zeros(N,N) for i in range(robo.NL)]
    
def init_vec(robo,N = 3,ext=0):
    """Generates a list of vectors. Size of the list is number of links.
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    N : int, optional
        size of the vectors, default is 3
    ext : int, optional
        additional vector instances over number of links
        
    Returns
    =======
    list of Matrices Nx1
    """
    return [zeros(N,1) for i in range(robo.NL+ext)]

def init_scalar(robo):
    """Generates a list of vars. Size of the list is number of links."""
    return [0 for i in range(robo.NL)]

def init_w(robo):
    """Generates a list of vectors for angular velocities.
    Size of the list is number of links + 1.
    The last vector is the base angular velocity"""
    w = init_vec(robo)
    w.append(robo.w0)
    return w
   
def init_wv_dot(robo):
    """Generates lists of vectors for angular and linear accelerations.
    Size of the list is number of links + 1.
    The last vector is the base angular velocity
            
    Returns
    =======
    vdot : list of Matrices 3x1
    wdot : list of Matrices 3x1
    """
    wdot = init_vec(robo)
    wdot.append(robo.wdot0)
    vdot = init_vec(robo)
    vdot.append(robo.vdot0 - robo.G)
    return wdot,vdot
    
def compute_transform(robo,symo,j,antRj,antPj):
    """Internal function. Computes rotation matrix and translation vector
    of ant_T_j homogenuous transform. Does the trigonometric subsctitution

    Notes
    =====
    antPj and antRj are the output parameters
    """
    robo.num[j] = robo.num[j]
    antTj = transform(robo,j)
    for angle,name in robo.get_angles(j):
        antTj = symo.trig_replace(antTj,angle,name)
    antRj[j] = symo.mat_replace(get_r(antTj),'A',robo.num[j])
    antPj[j] = symo.mat_replace(get_p(antTj),'L',robo.num[j])
    
def compute_screw_transform(robo,symo,j,antRj,antPj,jTant):
    """Internal function. Computes the screw transformation matrix
    between ant[j] and j frames.
    
    Notes
    =====
    jTant is an output parameter
    """
    jRant = antRj[j].T 
    ET = symo.mat_replace( - jRant*hat(antPj[j]),'JPR',robo.num[j])
    jTant[j] = (Matrix([jRant.row_join(ET),
                        zeros(3,3).row_join(jRant)]))
                        
def compute_twist(robo,symo,j,antRj,antPj,w,wdot,U,vdot):
    """Internal function. Computes angular velocity, auxiliary U matrix and 
    linear and angular accelerations.
            
    Notes
    =====
    w, wdot, U, vdot are the output parameters
    """
    jRant = antRj[j].T
    qdj = Matrix([0,0,robo.qdot[j]])
    qddj = Matrix([0,0,robo.qddot[j]])
    wi = symo.mat_replace(jRant*w[robo.ant[j]],'WI',robo.num[j])
    w[j] = symo.mat_replace(wi + (1 - robo.sigma[j])*qdj,'W',robo.num[j])
    wdot[j] = jRant*wdot[robo.ant[j]] + (1 - robo.sigma[j])*(qddj + hat(wi)*qdj)
    symo.mat_replace(wdot[j],'WP',robo.num[j])
    DV = product_combinations(w[j])
    symo.mat_replace(DV,'DV',robo.num[j]) 
    hatw_hatw = Matrix([[-DV[3]-DV[5],DV[1],DV[2]],
                        [DV[1],-DV[5]-DV[0],DV[4]],
                        [DV[2],DV[4],-DV[3]-DV[0]]])
    U[j] = hatw_hatw + hat(wdot[j])
    symo.mat_replace(U[j],'U',robo.num[j])    
    vsp = vdot[robo.ant[j]] + U[robo.ant[j]]*antPj[j]
    symo.mat_replace(vsp,'VSP',robo.num[j])
    vdot[j] = jRant*vsp + robo.sigma[j]*(qddj + 2*hat(wi)*qdj)
    symo.mat_replace(vdot[j],'VP',robo.num[j])  
                                
def compute_omega(robo,symo,j,antRj,w,wi): 
    """Internal function. Computes angular velocity of jth frame and
    projection of the antecedent frame's angular velocity
    
    Notes
    =====
    w, wi, U, vdot are the output parameters
    """
    
    jRant = antRj[j].T
    qdj = Matrix([0,0,robo.qdot[j]])
    wi[j] = symo.mat_replace(jRant*w[robo.ant[j]],'WI',robo.num[j])
    w[j] = symo.mat_replace(wi[j] + (1 - robo.sigma[j])*qdj,'W',robo.num[j])
                             
def compute_wrench(robo,symo,j,w,wdot,U,vdot,F,N):
    """Internal function. Computes total wrench (torques and forces)
    of the link j
    
    Notes
    =====
    F, N are the output parameters
    """
    F[j] = robo.M[j]*vdot[j] + U[j]*robo.MS[j]
    symo.mat_replace(F[j],'F',robo.num[j])
    Psi = symo.mat_replace(robo.J[j]*w[j],'PSI',robo.num[j])
    N[j] = robo.J[j]*wdot[j] + hat(w[j])*Psi
    symo.mat_replace(N[j],'No',robo.num[j])                        

def compute_joint_wrench(robo,symo,j,antRj,antPj,vdot,Fjnt,Njnt,F,N,Fex,Nex):
    """Internal function. Computes wrench (torques and forces)
    of the joint j
    
    Notes
    =====
    Fjnt, Njnt, Fex, Nex are the output parameters
    """
    Fjnt[j] = symo.mat_replace(F[j] + Fex[j],'E',robo.num[j])
    Njnt[j] = N[j] + Nex[j] + hat(robo.MS[j])*vdot[j]
    symo.mat_replace(Njnt[j],'N',robo.num[j])
    f_ant = symo.mat_replace(antRj[j]*Fjnt[j],'FDI',robo.num[j])
    if robo.ant[j] != - 1:
        Fex[robo.ant[j]] += f_ant
        Nex[robo.ant[j]] += antRj[j]*Njnt[j] + hat(antPj[j])*f_ant

def compute_torque(robo,symo,j,Fjnt,Njnt,name = 'GAM'):
    """Internal function. Computes actuation torques - projection of
    joint wrench on the joint axis
    """
    if robo.sigma[j] != 2:
        tau = (robo.sigma[j]*Fjnt[j] + (1 - robo.sigma[j])*Njnt[j])
        tau_total = tau[2] + robo.fric_s(j) + robo.fric_v(j) + robo.tau_ia(j)
        symo.replace(tau_total,name,robo.num[j],forced = True)

def inertia_spatial(J,MS,M):
    return Matrix([(M*eye(3)).row_join(hat(MS).T),hat(MS).row_join(J)])

def compute_joint_torque_deriv(symo,param,arg,index):
    """Internal function. Computes joint reactive torques
    in case if the parameter is 1
    
    Parameters
    ==========
    symo : Symoro
        symbol manager
    param : var
        Dynamic parameter
    arg : var
        The real torque is equal to arg*param
    index : strig
        identifies the parameter in the sybstituted symbol's name
    """
    if param != Integer(0) and arg != Integer(0):
        index = str(index) + str(param)
        symo.replace(arg,'DG',index,forced=True)

def compute_beta(robo,symo,j,w,beta_star):
    """Internal function. Computes link's wrench when 
    the joint accelerations are zero
    
    Notes
    =====
    beta_star is the output parameter
    """
    E1 = symo.mat_replace(robo.J[j]*w[j],'JW',robo.num[j])
    E2 = symo.mat_replace(hat(w[j])*E1,'KW',robo.num[j])
    E3 = hat(w[j])*robo.MS[j]
    E4 = symo.mat_replace(hat(w[j])*E3,'SW',robo.num[j])
    E5 = - robo.Nex[j] - E2
    E6 = - robo.Fex[j] - E4
    beta_star[j] = Matrix([E6,E5])

def compute_link_acc(robo,symo,j,antRj,antPj,link_acc,w,wi):
    """Internal function. Computes link's accelerations when 
    the joint accelerations are zero
    
    Notes
    =====
    link_acc is the output parameter
    """
    E1 = symo.mat_replace(hat(wi[j])*Matrix([0,0,robo.qdot[j]]),
                             'WQ',robo.num[j])
    E2 = (1 - robo.sigma[j])*E1
    E3 = 2*robo.sigma[j]*E1
    E4 = hat(w[robo.ant[j]])*antPj[j]
    E5 = hat(w[robo.ant[j]])*E4
    E6 = antRj[j].T*E5
    E7 = symo.mat_replace(E6 + E3,'LW',robo.num[j])
    link_acc[j] = Matrix([E7,E2])

def replace_beta_J_star(robo,symo,j,grandJ,beta_star):
    """Internal function. Makes symbol substitution in beta_star and grandJ
    """
    grandJ[j] = symo.mat_replace(grandJ[j],'MJE',robo.num[j],symmet = True)
    beta_star[j] = symo.mat_replace(beta_star[j],'VBE',robo.num[j])
  
def compute_Tau(robo,symo,j,grandJ,beta_star,jaj,juj,H_inv,Tau):
    """Internal function. Computes intermediat dynamic variables
    
    Notes
    =====
    H_inv and Tau are the output parameters
    """
    Jstar_jaj = grandJ[j]*jaj[j]
    if robo.sigma[j] == 2:
        Tau[j] = 0
    else:
        H_inv[j] = 1 / (jaj[j].dot(Jstar_jaj) + robo.IA[j])
        H_inv[j] = symo.replace(H_inv[j],'JD',robo.num[j])
        juj[j] = Jstar_jaj*H_inv[j]
        symo.mat_replace(juj[j],'JU',robo.num[j])
        joint_friction = robo.fric_s(j) + robo.fric_v(j)
        Tau[j] = jaj[j].dot(beta_star[j]) + robo.GAM[j] - joint_friction
        Tau[j] = symo.replace(Tau[j],'GW',robo.num[j])

def compute_beta_J_star(robo,symo,j,grandJ,jaj,juj,Tau,
                        beta_star,jTant,link_acc):
    """Internal function. Computes intermediat dynamic variables
    
    Notes
    =====
    grandJ and beta_star are the output parameters
    """
    Jstar_jaj = grandJ[j]*jaj[j]
    grandK = symo.mat_replace(grandJ[j] - juj[j]*Jstar_jaj.T,'GK',robo.num[j])
    E1 = symo.mat_replace(grandK*link_acc[j],'NG',robo.num[j])
    E3 = symo.mat_replace(E1 + Tau[j]*juj[j],'VS',robo.num[j])
    alpha = symo.mat_replace(E3 - beta_star[j],'AP',robo.num[j])
    E4 = symo.mat_replace(jTant[j].T*grandK,'GX',robo.num[j])
    E5 = symo.mat_replace(E4*jTant[j],'TKT',robo.num[j],symmet = True)
    grandJ[robo.ant[j]] += E5
    beta_star[robo.ant[j]] -= jTant[j].T*alpha
   
def compute_acceleration(robo,symo,j,jTant,grandVp,juj,H_inv,jaj,Tau,link_acc):
    """Internal function. Computes joint accelerations and links' twists
    
    Notes
    =====
    grandVp is the output parameter
    """
    grandR = symo.mat_replace(jTant[j]*grandVp[robo.ant[j]] + link_acc[j],
                              'VR',robo.num[j])       
    E1 = symo.replace(juj[j].dot(grandR),'GU',robo.num[j])
    if robo.sigma[j] == 2: 
        qddot = 0
    else:
        qddot = H_inv[j]*Tau[j] -  E1       
    qddot = symo.replace(qddot,"QDP",robo.num[j],forced = True)
    grandVp[j] = (grandR + qddot*jaj[j]) 
    grandVp[j][3:,0] = symo.mat_replace(grandVp[j][3:,0],'WP',robo.num[j])
    grandVp[j][:3,0] = symo.mat_replace(grandVp[j][:3,0],'VP',robo.num[j])
    
def compute_coupled_forces(robo,symo,j,grandVp,grandJ,beta_star):
    """Internal function.
    """
    E3 = symo.mat_replace(grandJ[j]*grandVp[j],'DY',robo.num[j])
    couplforce = E3 - beta_star[j];
    symo.mat_replace(couplforce[3:,0],'N',robo.num[j])
    symo.mat_replace(couplforce[:3,0],'E',robo.num[j])
    
def replace_Jplus(robo,symo,j,Jplus,MSplus,Mplus):
    """Internal function. Makes symbol substitutions inertia parameters
    """
    symo.mat_replace(Jplus[j],'JP',robo.num[j])
    symo.mat_replace(MSplus[j],'MSP',robo.num[j])
    Mplus[j] = symo.replace(Mplus[j],'MP',robo.num[j])

def compute_Jplus(robo,symo,j,antRj,antPj,Jplus,MSplus,Mplus,AJE1):
    """Internal function. Computes inertia parameters of composed link
    
    Notes
    =====
    Jplus, MSplus, Mplus are the output parameters
    """
    hat_antPj = hat(antPj[j])
    antMSj = symo.mat_replace(antRj[j]*MSplus[j],'AS',robo.num[j])
    E1 = symo.mat_replace(antRj[j]*Jplus[j],'AJ',robo.num[j])
    AJE1[j] = E1[:,2]
    E2 = symo.mat_replace(E1*antRj[j].T,'AJA',robo.num[j])
    E3 = symo.mat_replace(hat_antPj*hat(antMSj),'PAS',robo.num[j])    
    Jplus[robo.ant[j]] += E2 - (E3 + E3.T) + hat_antPj*hat_antPj.T*Mplus[j]
    MSplus[robo.ant[j]] += antMSj + antPj[j]*Mplus[j]
    Mplus[robo.ant[j]] += Mplus[j]

def compute_A_diagonal(robo,symo,j,Jplus,MSplus,Mplus,f,n,A):
    """Internal function. Computes diagonal elements
    of the inertia matrix
    
    Notes
    =====
    f, n, A are the output parameters
    """
    if robo.sigma[j]==0:
        f[j]=Matrix([ - MSplus[j][1],MSplus[j][0],0])
        n[j]=Jplus[j][:,2]                       
        A[j,j]=Jplus[j][2,2] + robo.IA[j]        
    elif robo.sigma[j] == 1:
        f[j]=Matrix([0,0,Mplus[j]])
        n[j]=Matrix([MSplus[j][1], - MSplus[j][0],0])
        A[j,j]=Mplus[j] + robo.IA[j]  
    symo.mat_replace(f[j],'E' + chars[j],robo.num[j])
    symo.mat_replace(n[j],'N' + chars[j],robo.num[j])
     
def compute_A_triangle(robo,symo,j,k,ka,antRj,antPj,f,n,A,AJE1):
    """Internal function. Computes elements below and above diagonal
    of the inertia matrix
    
    Notes
    =====
    f, n, A are the output parameters
    """
    f[ka] = antRj[k]*f[k]
    if k == j and robo.sigma[j] ==0:
        n[ka] = AJE1[k] + hat(antPj[k])*f[k]
    else:
        n[ka] = antRj[k]*n[k] + hat(antPj[k])*f[k]
    if ka == - 1:
        symo.mat_replace(f[ka],'AV0')
        symo.mat_replace(n[ka],'AW0')
    else:
        symo.mat_replace(f[ka],'E' + chars[j],robo.num[ka])
        symo.mat_replace(n[ka],'N' + chars[j],robo.num[ka]) 
        if robo.sigma[ka] == 0:
            A[j,ka] = n[ka][2]
        elif robo.sigma[ka] == 1:
            A[j,ka] = f[ka][2]
        A[ka,j] = A[j,ka]

#TODO:Finish base parameters computation
def base_paremeters(robo_orig):
    """Computes grouped inertia parameters. New parametrization contains
    less parameters but generates the same dynamic model
    
    Parameters
    ==========
    robo : Robot
        Instance of robot description container
    
    Returns
    =======
    symo.sydi : dictionary
        Dictionary with the information of all the sybstitution
    """
    robo = copy(robo_orig)
    lam = [0 for i in range(robo.NL)]
    antRj = init_mat(robo)
    antPj = init_vec(robo)
    symo = Symoro()
    symo.file_open(robo,'regp')
    symo.write_dynam_param(robo,'Base parameters computation')
    for j in range(robo.NL):
        compute_transform(robo,symo,j,antRj,antPj)
    for j in reversed(range(robo.NL)):
        replace_dynam_param(robo,symo,j)
        if robo.ant[j] == - 1:
            pass
        elif robo.sigma[j] == 0:
            #general grouping
            compute_lambda(robo,symo,j,antRj,antPj,lam)
            group_param_rot(robo,symo,j,lam)
            #special grouping
            pass
        elif robo.sigma[j] == 1:            
            #general grouping
            group_param_prism(robo,symo,j,antRj)
            #special grouping
            pass
        elif robo.sigma[j] == 2:
            compute_lambda(robo,symo,j,antRj,antPj)
            group_param_fix(robo,symo,j,lam)
            pass       
        pass
    symo.write_line()
    symo.write_inert_param(robo,robo.name + ' grouped inertia parameters.')
    symo.file_out.close()
    return robo,symo.sydi
        
def vec_mut_J(v,u):
    """Internal function. Needed for inertia parameters transformation
    
    Parameters
    ==========
    v, u : Matrix 3x1
        two axis vectors
    Returns : Matrix 6x1
    """
    return Matrix([v[0]*u[0],v[0]*u[1],v[0]*u[2],v[1]*u[1],v[1]*u[2],v[2]*u[2]])
    
def vec_mut_MS(v,P):
    """Internal function. Needed for inertia parameters transformation
    
    Parameters
    ==========
    v : Matrix 3x1
        axis vector
    P : Matrix 3x1
        position vector
    
    Returns : Matrix 6x1
    """
    U = - hat(v)*hat(P)
    return Matrix([2*U[0,0],U[0,1] + U[1,0],U[0,2] + U[2,0],\
                2*U[1,1],U[1,2] + U[2,1],2*U[2,2]])

def vec_mut_M(P):
    """Internal function. Needed for inertia parameters transformation
    
    Parameters
    ==========
    P : Matrix 3x1
        position vector
    
    Returns : Matrix 6x1
    """
    U = - hat(P)*hat(P)
    return  Matrix([U[0,0],U[0,1],U[0,2],U[1,1],U[1,2],U[2,2]])   
    
def compute_lambda(robo,symo,j,antRj,antPj,lam):
    """Internal function. Computes the inertia parameters transformation matrix
    
    Notes
    =====
    lam is the output paramete
    """
    lamJJ_list = []
    lamJMS_list = []
    for e1 in range(3):
        for e2 in range(e1,3):
            u = vec_mut_J(antRj[j][:,e1],antRj[j][:,e2])
            if e1 != e2:
                u += vec_mut_J(antRj[j][:,e2],antRj[j][:,e1])
            lamJJ_list.append(u.T)
    for e1 in range(3):
        v = vec_mut_MS(antRj[j][:,e1],antPj[j])
        lamJMS_list.append(v.T)
    lamJJ = Matrix(lamJJ_list).T #,'LamJ',robo.num[j])
    lamJMS = symo.mat_replace(Matrix(lamJMS_list).T,'LamMS',robo.num[j])
    lamJM = symo.mat_replace(vec_mut_M(antPj[j]),'LamM',robo.num[j])
    lamJ = lamJJ.row_join(lamJMS).row_join(lamJM)
    lamMS = zeros(3,6).row_join(antRj[j]).row_join(antPj[j])
    lamM = zeros(1,10)
    lamM[9] = 1
    lam[j] = Matrix([lamJ,lamMS,lamM])

def group_param_rot(robo,symo,j,lam):
    """Internal function. Groups inertia parameters according to the
    general rule for a rotational joint.
    
    Notes
    =====
    robo is the output paramete
    """
    Kj = robo.get_dynam_param(j)
    Kant = robo.get_dynam_param(robo.ant[j])
    lam03 = lam[j][:,0] + lam[j][:,3]
    symo.mat_C2S2_simp(lam03,robo.num[j])
    Kant += lam03*Kj[3] + lam[j][:,8]*Kj[8] + lam[j][:,9]*Kj[9]
    Kj[3] = 0   #YY
    Kj[8] = 0   #MZ  
    Kj[9] = 0   #M
    robo.put_dynam_param(Kant,robo.ant[j])
    robo.put_dynam_param(Kj,j)

def group_param_fix(robo,symo,j,lam):
    """Internal function. Groups inertia parameters according to the
    general rule for a fixed joint joint.
    
    Notes
    =====
    robo is the output paramete
    """
    Kj = robo.get_dynam_param(j)
    Kant = robo.get_dynam_param(robo.ant[j])
    Kant += lam[j]*Kj
    robo.put_dynam_param(Kant,robo.ant[j])
    robo.put_dynam_param(zeros(10,1),j)
    
def replace_dynam_param(robo,symo,j):
    """Internal function. Makes symbol substitutions for inertia parameters
    
    Notes
    =====
    robo is the output parameter
    """
    symo.mat_replace(robo.J[j],'JR',robo.num[j])
    symo.mat_replace(robo.MS[j],'MSR',robo.num[j])
    robo.M[j] = symo.replace(robo.M[j],'MR',robo.num[j])
    
def group_param_prism(robo,symo,j,antRj):
    """Internal function. Groups inertia parameters according to the
    general rule for a prismatic joint.
    
    Notes
    =====
    robo is the output paramete
    """
    antJj = antRj[j]*robo.J[j]*antRj[j].T
    robo.J[robo.ant[j]] += antJj
    robo.J[j] = zeros(3,3)  

robo = RX90()
     
print 'Inverse dynamic model using Newton - Euler Algorith'
inverse_dynamic_NE(robo)

print 'Pseudo forces using Newton - Euler Algorith'
pseudo_force_NE(robo)

print 'Direct dynamic model using Newton - Euler Algorith'
direct_dynamic_NE(robo)

print 'Dynamic identification model using Newton - Euler Algorith'
dynamic_identification_NE(robo)

print 'Inertia Matrix using composite links'
inertia_matrix(robo)

print 'Base parameters computation'
base_paremeters(robo)
