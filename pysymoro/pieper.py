# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module of SYMORO package provides symbolic
solutions for inverse geometric problem using Pieper Method.
"""


from sympy import var, sin, cos, eye, atan2, sqrt, pi
from sympy import Matrix, Symbol, Expr, trigsimp, zeros, ones
from numpy import array, dot

from symoroutils import symbolmgr, tools
from pysymoro import invdata
from pysymoro.invdata import solve_position
from pysymoro.invdata import solve_orientation
from pysymoro.invdata import solve_orientation_prismatic
from pysymoro.invdata import solve_position_prismatic


# dictionary for joint combinations out of the pieper ones
# coding: (i,j,k):type
joint_com = {
    (0, 0, 0): 0, (0, 0, 1): 1, (0, 1, 0): 2,
    (0, 1, 1): 3, (1, 0, 0): 4, (1, 0, 1): 5,
    (1, 1, 0): 6, (1, 1, 1): 7
}


def _pieper_solve(robo, symo):
    """
    Main function of the Pieper solutions.

    Parameters:
    ===========
    1) robo: Parameter that gives us access to the parameters of the robot.
    2) symo: Instance that gives us access for the symbolic handling.
    """
    if robo.structure == tools.SIMPLE:
        [pieper_joints, bools] = _look_for_case_simple(robo,symo)
        [bool_fail,bool_prism,bool_spherical] = bools
        if bool_fail == 0:
            [com_key, X_joints] = _X_joints(robo,symo,pieper_joints)
        else:
            return
        if bool_spherical == 1:
            m = pieper_joints[1]
            if m == 2:
                RRRXXX(robo,symo,com_key,X_joints,pieper_joints)
            elif m == 3:
                XRRRXX(robo,symo,com_key,X_joints,pieper_joints)
            elif m == 4:
                XXRRRX(robo,symo,com_key,X_joints,pieper_joints)
            elif m == 5:
                XXXRRR(robo,symo,com_key,X_joints,pieper_joints)
        elif bool_prism == 1:
            Prismatic(robo,symo,X_joints,pieper_joints)
    elif (robo.structure == tools.TREE) or (robo.structure == tools.CLOSED_LOOP):
        [bools, pieper_branches, pieper_joints_all, X_joints_all, com_key] = _look_for_case_tree(robo,symo)
        [bool_fail, bool_prism, bool_spherical] = bools
        if sum(bool_fail) == len(bool_fail):
            return
        X_joints = [X_joints_all[x:x+3] for x in range(0, len(X_joints_all), 3)]
        pieper_joints = [pieper_joints_all[x:x+3] for x in range(0, len(pieper_joints_all), 3)]
        for i in range(len(pieper_branches)):
            symo.write_line("# Solution for branch {0} :".format(i+1))
            symo.write_line("#============================= \r\n\r\n")
            for j in range(len(bool_spherical)):
                if bool_spherical[j] == 1:
                    m = pieper_joints[i][1]
                    bool_spherical[j] = 0
                    if m == 2:
                        RRRXXX(robo,symo,com_key[i],X_joints[i],pieper_joints[i])
                        break
                    elif m == 3:
                        XRRRXX(robo,symo,com_key[i],X_joints[i],pieper_joints[i])
                        break
                    elif m == 4:
                        XXRRRX(robo,symo,com_key[i],X_joints[i],pieper_joints[i])
                        break
                    elif m == 5:
                        XXXRRR(robo,symo,com_key[i],X_joints[i],pieper_joints[i])
                        break
                elif bool_prism[j] == 1:
                    bool_prism[j] = 0
                    Prismatic(robo,symo,X_joints[i],pieper_joints[i])
                    break
    return


def _look_for_case_simple(robo, symo):
    """
    Function that locates if a serial robot can be solved by PIEPER METHOD.

    Parameters:
    ===========
    """
    try_paul_str = ("\r\n\r\n# This robot cannot be ") + \
        ("solved by PIEPER METHOD. Try Paul method. \r\n\r\n")
    symo.write_line("# Pieper Method data for the given robot: ")
    [bool_fail, bool_prism, bool_spherical] = [0,0,0]
    pieper_joints = []
    if robo.nj > 6:
        symo.write_line(
            ("\r\n\r\n# This robot cannot be solved by PIEPER METHOD.") + \
            (" The robot is redundant. \r\n\r\n")
        )
        bool_fail = 1
    elif robo.nj < 6:
        symo.write_line(try_paul_str)
        symo.write_line("# --The robot has less than 6-DOF.--")
        bool_fail = 1
    else:
        if sum(robo.sigma[1:len(robo.sigma)]) > 3:
            symo.write_line(
                ("\r\n\r\n# This robot cannot be solved by ") + \
                ("PIEPER METHOD, because it is redundant. ") + \
                ("(more than three prismatic joints)  \r\n\r\n")
            )
            bool_fail = 1
        elif sum(robo.sigma[1:len(robo.sigma)]) == 3:
            bool_prism = 1
            for joint in range(1, len(robo.sigma)):
                if robo.sigma[joint] == 1:
                    pieper_joints.append(joint)
            symo.write_line(
                ("\r\n\r\n# PIEPER METHOD: Decoupled 6-DOF robot ") + \
                ("with 3 prismatic joints positioned at:") + \
                ("{0}".format(pieper_joints)) + \
                ("\r\n\r\n")
            )
        else:
            for m in range(2, len(robo.sigma)-1):
                if (robo.sigma[m-1] == 0) and (robo.sigma[m] == 0) and (robo.sigma[m+1] == 0):
                    if (robo.d[m] == 0) and (robo.d[m+1] == 0) and (robo.r[m] == 0):
                        if (sin(robo.alpha[m]) != 0) and (sin(robo.alpha[m+1]) != 0):
                            pieper_joints = [m-1,m,m+1]
                            symo.write_line(
                                ("\r\n\r\n# PIEPER METHOD: Decoupled ") + \
                                ("6-DOF robot with a spherical joint ") + \
                                ("composed by joints: ") + \
                                ("{0}".format(pieper_joints)) + \
                                ("\r\n\r\n")
                            )
                            bool_spherical = 1
                            break
                        elif m == 5:
                            bool_fail = 1
                            symo.write_line(try_paul_str)
                    elif m == 5:
                        bool_fail = 1
                        symo.write_line(try_paul_str)
                elif m == 5:
                    bool_fail = 1
                    symo.write_line(try_paul_str)
    bools = [bool_fail, bool_prism, bool_spherical]
    return pieper_joints, bools


def _look_for_case_tree(robo,symo):
    try_paul_str =  ("\r\n\r\n# Branch {0}") + \
        ("of the robot cannot be solved by PIEPER METHOD. ") + \
        ("This branch has less than 6 joints. ") + \
        ("Try Paul Method \r\n\r\n")
    End_joints = []
    X_joints = []
    pieper_joints = []
    j = []
    for i in range(robo.NJ+1):
        j.append(i)
    for joint in range(1, len(robo.ant)):
        if j[joint] not in robo.ant and robo.sigma[joint] != 2:
            End_joints.append(j[joint])
            branches = len(End_joints)
    symo.write_line("# The tree structure robot has {0} branches.".format(branches))
    [bool_fail, bool_prism, bool_spherical] = [zeros(1, branches), zeros(1, branches), zeros(1, branches)]
    [pieper_branches, com_key] = [999*ones(1, branches), 999*ones(1, branches)]
    num_prism = zeros(1, branches)
    for i in range(branches):
        f_joint = End_joints[i]
        c = 0
        globals()["bran"+str(i)] = [0]*max(robo.ant)
        while f_joint != 0:
            globals()["bran"+str(i)][c] = f_joint
            f_joint = robo.ant[f_joint]
            c += 1
        globals()["bran"+str(i)] = [x for x in globals()["bran"+str(i)] if x != 0]
        globals()["bran"+str(i)] = globals()["bran"+str(i)][::-1]
        if len(globals()["bran"+str(i)]) > 6:
            symo.write_line(
                ("\r\n\r\n# Branch {0} ".format(i+1)) + \
                ("of the robot cannot be solved by PIEPER METHOD. ") + \
                ("This branch is redundant. \r\n\r\n")
            )
            bool_fail[i] = 1
        elif len(globals()["bran"+str(i)]) < 6:
            symo.write_line(try_paul_str.format(i+1))
            bool_fail[i] = 1
        else:
            bool_fail[i] = 0
            for k in range(len(globals()["bran"+str(i)])):
                num_prism[i] = num_prism[i] + robo.sigma[globals()["bran"+str(i)][k]]
            if num_prism[i] > 3:
                symo.write_line(
                    ("\r\n\r\n# Branch {0} ".format(i+1)) + \
                    ("of the robot cannot be solved by PIEPER METHOD. ") + \
                    ("It is redundant (more than 3 prismatic joints). \r\n\r\n")
                )
                bool_fail[i] = 1
            elif num_prism[i] == 3:
                pieper_branches[i] = i
                bool_prism[i] = 1
                p_joints = []
                joints = []
                for joint in range(len(globals()["bran"+str(i)])):
                    if robo.sigma[globals()["bran"+str(i)][joint]] == 1:
                        pieper_joints.append(globals()["bran"+str(i)][joint])
                        p_joints.append(globals()["bran"+str(i)][joint])
                    else:
                        X_joints.append(globals()["bran"+str(i)][joint])
                        joints.append(robo.sigma[globals()["bran"+str(i)][joint]])
                        joint_type = tuple(joints)
                    if joint_type in joint_com:
                        com_key[i] = joint_com[joint_type]
                symo.write_line(
                    ("\r\n\r\n# PIEPER METHOD: Branch {0}".format(i+1)) + \
                    (" is decoupled with 3 prismatic joints ") + \
                    ("positioned at {0} \r\n\r\n".format(p_joints))
                )
            else:
                for m in range(2, len(globals()["bran"+str(i)])):
                    if (robo.sigma[globals()["bran"+str(i)][m-1]] == 0) \
                        and (robo.sigma[globals()["bran"+str(i)][m]] == 0) \
                        and (robo.sigma[globals()["bran"+str(i)][m+1]] == 0):
                        if (robo.d[globals()["bran"+str(i)][m]] == 0) \
                            and (robo.d[globals()["bran"+str(i)][m+1]] == 0) \
                            and (robo.r[globals()["bran"+str(i)][m]] == 0):
                            if (sin(robo.alpha[globals()["bran"+str(i)][m]]) != 0) \
                                and (sin(robo.alpha[globals()["bran"+str(i)][m+1]]) != 0):
                                pieper_branches[i] = i
                                pieper_joints = [
                                    globals()["bran"+str(i)][m-1],
                                    globals()["bran"+str(i)][m],
                                    globals()["bran"+str(i)][m+1]
                                ]
                                joints = []
                                joint = globals()["bran"+str(i)]
                                for ji in range(len(joint)):
                                    if joint[ji] not in pieper_joints:
                                        X_joints.append(globals()["bran"+str(i)][ji])
                                        joints.append(robo.sigma[globals()["bran"+str(i)][ji]])
                                        joint_type = tuple(joints)
                                    if joint_type in joint_com:
                                        com_key[i] = joint_com[joint_type]
                                symo.write_line(
                                    ("\r\n\r\n# PIEPER METHOD: ") + \
                                    ("Branch{0}".format(i+1)) + \
                                    (" is decoupled with a ") + \
                                    ("spherical joint composed by ") + \
                                    ("joints {0} \r\n\r\n".format(pieper_joints))
                                )
                                bool_spherical[i] = 1
                                break
                            elif m == 5:
                                bool_fail[i] = 1
                                symo.write_line(try_paul_str.format(i+1))
                        elif m == 5:
                            bool_fail[i] = 1
                            symo.write_line(try_paul_str.format(i+1))
                    elif m == 5:
                        bool_fail[i] = 1
                        symo.write_line(try_paul_str.format(i+1))
    bools = [bool_fail, bool_prism, bool_spherical]
    pieper_branches = [x for x in pieper_branches if x != 999]
    com_key = [x for x in com_key if x != 999]
    return bools, pieper_branches, pieper_joints, X_joints, com_key


def igm_pieper(robo, T_ref, n):
    symo = symbolmgr.SymbolManager()
    symo.file_open(robo, 'pieper')
    symo.write_params_table(robo, 'Inverse Geometrix Model for frame %s' % n)
    _pieper_solve(robo, symo)
    symo.file_close()
    return symo


def _X_joints(robo, symo, pieper_joints):
    """
    Function that takes the already identified pieper_joints and returns the other X_joints.

    Parameters:
    ===========
    """
    # pieper_joints will be a vector containing the position of the
    # joints in the parameters table [pieper_joints = (q1,q2,q3)]
    joints = []                         # Empty vector to append the type of the joints
    X_joints = []
    for joint in range(1, len(robo.sigma)):
        if joint not in pieper_joints:
            X_joints.append(joint)
            joints.append(robo.sigma[joint])
            joint_type = tuple(joints)        # Convert from list to tuple in order to use the vector
    # Identify the combination of the X joints
    if joint_type in joint_com: com_key = joint_com[joint_type]

    return com_key, X_joints


def XXXRRR(robo, symo, com_key, X_joints, pieper_joints):
    """
    Function that prints the symbolic solution of this decoupled robot case using PIEPER METHOD.

    Parameters:
    ===========
    1) X_joints: Joints that are out of spherical wrist
    2) pieper_joints: Joints that the spherical wrist is composed of.
    """
    ## P36_1 = d4
    P36_1 = robo.d[4]
    ## P36_2 = -r4*Sa4
    P36_2 = -sin(robo.alpha[4])*robo.r[4]
    ## P36_3 = r4*Ca4
    P36_3 = cos(robo.alpha[4])*robo.r[4]
    [px,py,pz] = invdata.T_GENERAL[:3,3]
    g = Matrix([px,py,pz,1])
    ## fc = [d4; -Ca3*Sa4*r4; -Sa3*Sa4*r4]
    fc = Matrix([P36_1, cos(robo.alpha[3])*P36_2, sin(robo.alpha[3])*P36_2, 0])
    ## fs = [Sa4*r4; Ca3*d4; Sa3*d4]
    fs = Matrix([-P36_2, cos(robo.alpha[3])*P36_1, sin(robo.alpha[3])*P36_1, 0])
    ## fr = [0; -Sa3; Ca3]
    fr = Matrix([0, -sin(robo.alpha[3]), cos(robo.alpha[3]), 0])
    ## f0 = [d3; -Sa3*Ca4*r4; Ca3*Ca4*r4]
    f0 = Matrix([robo.d[3], -sin(robo.alpha[3])*P36_3, cos(robo.alpha[3])*P36_3, 1])
    # Position Equations First
    solve_position(robo, symo, com_key, X_joints, fc, fs, fr, f0, g)
    # Then Orientation Equations
    solve_orientation(robo, symo, pieper_joints)
    return


def XRRRXX(robo, symo, com_key, X_joints, pieper_joints):
    """
    Function that prints the symbolic solution of this decoupled robot case using PIEPER METHOD.

    Parameters:
    ===========
    1) X_joints: Joints that are out of spherical wrist
    2) pieper_joints: Joints that the spherical wrist is composed of.
    """
    g = Matrix([robo.d[2], -sin(robo.alpha[2])*robo.r[2], cos(robo.alpha[2])*robo.r[2], 1])
    A0 = array(invdata.T_GENERAL[:3, :3])
    P0 = array(invdata.T_GENERAL[:3, 3])
    ## fc = A0*[d5-d6; r5*(Ca5*Sa6-Ca6*Sa5); 0]
    fc = dot(
        A0,
        [robo.d[5]-robo.d[6],
        robo.r[5]*(cos(robo.alpha[5])*sin(robo.alpha[6])-cos(robo.alpha[6])*sin(robo.alpha[5])),
        0]
    )
    fc = Matrix([fc[0], fc[1], fc[2], 0])
    ## fs = A0*[r5*(Ca5*Sa6-Ca6*Sa5); d6-d5; 0]
    fs = dot(
        A0,
        [robo.r[5]*(cos(robo.alpha[5])*sin(robo.alpha[6])-cos(robo.alpha[6])*sin(robo.alpha[5])),
        robo.d[6]-robo.d[5], 0]
    )
    fs = Matrix([fs[0], fs[1], fs[2], 0])
    ## fr = A0*[0; 0; -1]
    fr = dot(A0, [0, 0, -1])
    fr = Matrix([fr[0], fr[1], fr[2], 0])
    ## f0 = A0*[0; 0; r5*C(a5-a6)] + P0
    f0 = sum(dot(A0, [0, 0, robo.r[5]*cos(robo.alpha[5]-robo.alpha[6])]), P0)
    f0 = Matrix([f0[0][0], f0[1][0], f0[2][0], 1])
    # Position Equations First
    solve_position(robo, symo, com_key, X_joints, fc, fs, fr, f0, g)
    # Then Orientation Equations
    solve_orientation(robo, symo, pieper_joints)
    return


def XXRRRX(robo, symo, com_key, X_joints, pieper_joints):
    """
    Function that prints the symbolic solution of this decoupled robot case using PIEPER METHOD.

    Parameters:
    ===========
    1) X_joints: Joints that are out of spherical wrist
    2) pieper_joints: Joints that the spherical wrist is composed of.
    """
    g = Matrix([robo.d[3], -sin(robo.alpha[3])*robo.r[3], cos(robo.alpha[3])*robo.r[3], 1])
    A0 = array(invdata.T_GENERAL[:3, :3])
    P0 = array(invdata.T_GENERAL[:3, 3])
    fc = dot(A0, [-robo.d[6], -robo.r[5]*sin(robo.alpha[6]), 0])
    ## fc = A0*[-d6; -Sa6*r5; 0]
    fc = Matrix([fc[0], fc[1], fc[2], 0])
    symo.write_line(fc)
    ## fs = A0*[-Sa6*r5; d6; 0]
    fs = dot(A0, [-robo.r[5]*sin(robo.alpha[6]), robo.d[6], 0])
    fs = Matrix([fs[0], fs[1], fs[2], 0])
    ## fr = -A0*[0; 0; 1]
    fr = -dot(A0, [0, 0, 1])
    fr = Matrix([fr[0], fr[1], fr[2], 0])
    ## f0 = A0*[0; 0; -Ca6*r5] + P0
    f0 = sum(dot(A0, [0, 0, -robo.r[5]*cos(robo.alpha[6])]), P0)
    f0 = Matrix([f0[0][0], f0[1][0], f0[2][0], 1])
    # Position Equations First
    solve_position(robo, symo, com_key, X_joints, fc, fs, fr, f0, g)
    # Then Orientation Equations
    solve_orientation(robo, symo, pieper_joints)
    return


def RRRXXX(robo, symo, com_key, X_joints, pieper_joints):
    """
    Function that prints the symbolic solution of this decoupled robot case using PIEPER METHOD.

    Parameters:
    ===========
    1) X_joints: Joints that are out of spherical wrist
    2) pieper_joints: Joints that the spherical wrist is composed of.
    """
    ## P63_1 = -d4
    P63_1 = -robo.d[4]
    ## P63_2 = -r3*Sa4
    P63_2 = -sin(robo.alpha[4])*robo.r[3]
    ## P63_3 = -r3*Ca4
    P63_3 = -cos(robo.alpha[4])*robo.r[3]
    A0 = array(invdata.T_GENERAL[:3, :3])
    P0 = array(invdata.T_GENERAL[:3, 3])
    g = -dot(A0, P0)
    g = Matrix([g[0][0], g[1][0], g[2][0], 1])
    ## fc = [-d4; -Ca5*Sa4*r3; Sa5*Sa4*r3]
    fc = Matrix([P63_1, cos(robo.alpha[5])*P63_2, -sin(robo.alpha[5])*P63_2, 0])
    ## fs = [Sa4*r3; -Ca5*d4; Sa5*d4]
    fs = Matrix([-P63_2, cos(robo.alpha[5])*P63_1, -sin(robo.alpha[5])*P63_1, 0])
    ## fr = [0; Sa5; Ca5]
    fr = Matrix([0, sin(robo.alpha[5]), cos(robo.alpha[5]), 0])
    ## f0 = [-d5; -Sa5*Ca4*r3; -Ca5*Ca4*r3]
    f0 = Matrix([-robo.d[5], sin(robo.alpha[5])*P63_3, cos(robo.alpha[5])*P63_3, 1])
    # Position Equations First
    solve_position(robo, symo, com_key, X_joints, fc, fs, fr, f0, g)
    # Then Orientation Equations
    solve_orientation(robo, symo, pieper_joints)
    return


def Prismatic(robo, symo, X_joints, pieper_joints):
    """
    Function that prints the symbolic solution of this decoupled robot case using PIEPER METHOD.

    Parameters:
    ===========
    1) X_joints: The three revolute joints.
    2) pieper_joints: The three prismatic joints.
    """
    # Orientation Equations First
    solve_orientation_prismatic(robo, symo, X_joints)
    # Then Position Equations
    solve_position_prismatic(robo, symo, pieper_joints)
    return


