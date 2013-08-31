"""
This module of SYMORO package provides symbolic
solutions for inverse geompetric problem.

The core symbolic library is sympy.
Needed modules : symoro.py, geometry.py

ECN - ARIA1 2013
"""

from sympy import Matrix, var, Symbol, sin, cos, eye, atan2, sqrt, pi
from symoro import Symoro, Robot, ZERO, ONE, get_max_coef
from geometry import dgm
from heapq import heapify, heappop

EMPTY = var("EMPTY")

#dictionary for equation type classification
eq_dict = {(1, 0, 0):0, (0, 1, 0):1, (1, 1, 0):2,
                  (0, 2, 0):3, (0, 2, 1):4}
name_dict = ("type 1", "type 2, 3", "type 4, 5", "type 6, 7", "type 8")

#def print_eq(symo, A, B):
#    if isinstance(B, Expr):
#        B = symo.CS12_simp(B)
#    print "    " + str(A) + " = " + str(B)
#    return B

def paul_solve(robo, symo, nTm, n, m, known = set()):
    chain = robo.loop_chain(m,n)
    iTn = dgm(robo, symo, m, n, key = 'left', trig_subs = False)
    iTm = dgm(robo, symo, n, m, key = 'left', trig_subs = False)
    th_all = set(robo.theta[i] for i in chain if i >= 0 and robo.sigma[i] == 0)
    r_all = set(robo.r[i] for i in chain if i >= 0 and  robo.sigma[i] == 1)
    while True:
        repeat_search = False
        for i in reversed(chain):
            M_eq = iTn[(i, n)]*nTm - iTm[(i, m)]
            while True:
                cont_search = False
                eq_candidates = [list() for list_index in xrange(5)]
                for e1 in xrange(3):
                    for e2 in xrange(4):
                        if M_eq[e1, e2].has(EMPTY):
                            continue
                        eq = symo.unknown_sep(M_eq[e1, e2],known)#.expand()
                        th_vars = (eq.atoms(Symbol) & th_all) - known
                        if th_vars:
                            arg_sum = max(at.count_ops()-1 for at in eq.atoms(sin, cos)
                                        if not at.atoms(Symbol) & known)
                        else:
                            arg_sum = 0
                        r_vars = (eq.atoms(Symbol) & r_all) - known
                        eq_features = (len(r_vars), len(th_vars), arg_sum)
                        if eq_features in eq_dict:
                            eq_key = eq_dict[eq_features]
                            eq_pack = (eq, list(r_vars), list(th_vars))
                            eq_candidates[eq_key].append(eq_pack)
                cont_search |= try_solve_0(symo, eq_candidates[0], known)
                cont_search |= try_solve_1(symo, eq_candidates[1], known)
                cont_search |= try_solve_2(symo, eq_candidates[2] + eq_candidates[1], known)
                cont_search |= try_solve_3(symo, eq_candidates[3], known)
                cont_search |= try_solve_4(symo, eq_candidates[4], known)
                repeat_search |= cont_search
                if not cont_search or th_all | r_all <= known:
                    break
            if th_all | r_all <= known:
                break
        if not repeat_search or th_all | r_all <= known:
            break
    return known

def loop_solve(robo, symo, knowns = None):
    #TODO: rewrite
    q_vec = robo.get_q_vec()
    loops = []
    if knowns == None:
        knowns = set(q for i, q in enumerate(q_vec) if robo.mu[i] == 1)
    for i, j in robo.get_loop_terminals():
        chain = robo.loop_chain(i,j)
        knowns_ij = set(q_vec[i] for i in chain if q_vec[i] in knowns)
        unknowns_ij = set(q_vec[i] for i in chain if q_vec[i] not in knowns)
        loops.append([len(unknowns_ij),i,j,knowns_ij,unknowns_ij])
    while loops:
        heapify(loops)
        loop = heappop(loops)
        res_knowns = paul_solve(robo, symo, eye(4), *loop[1:4])
        for l in loops:
            found = l[4] & res_knowns
            l[3] |= found
            l[4] -= found
            l[0] = len(l[4])

def igm_Paul(robo, symo, T_ref):
    terminal = set(range(robo.NL))-set(robo.ant)
    for n in terminal:   #terminal:
        print "#solution for the frame", robo.num[n]
        paul_solve(robo, symo, T_ref, -1, n)

def try_solve_0(symo, eq_sys, known):
    res = False
    for eq, [r], th_names in eq_sys:
        X = get_max_coef(eq, r)
        if X != 0:
            Y = X*r - eq
            print "type 1"
            X = symo.replace(symo.CS12_simp(X), 'X', r)
            Y = symo.replace(symo.CS12_simp(Y), 'Y', r)
            symo.add_to_dict(r, Y/X)
            known.add(r)
            res = True
    return res

def try_solve_1(symo, eq_sys, known):
    res = False
    for i in xrange(len(eq_sys)):
        eqi, r_vari, [th_i ] = eq_sys[i]
        if th_i  in known:
            continue
        Xi, Yi, Zi, i_ok = get_coefs(eqi, sin(th_i ), cos(th_i ))
        i_ok &= sum([Xi == ZERO, Yi == ZERO, Zi == ZERO]) <= 1
        if not i_ok:
            continue
        j_ok = False
        for j in xrange(i+1, len(eq_sys)):
            eqj, r_varj, [th_j] = eq_sys[j]
            if th_i  == th_j:
                Xj, Yj, Zj, j_ok = get_coefs(eqj, sin(th_j), cos(th_j))
                j_ok &= (Xi*Yj != Xj*Yi)
                if j_ok:
                    break
        if j_ok:
            print "#Solving type 3"
            solve_type_3(symo, Xi, Yi, -Zi, Xj, Yj, -Zj, th_i)
        else:
            print "#Solving type 2"
            solve_type_2(symo, Xi, Yi, -Zi, th_i)
        known.add(th_i )
        res = True
    return res

def try_solve_2(symo, eq_sys, known):
    if all(len(r_var) == 0 for eq, r_var, ths in eq_sys):
        return False
    for i in xrange(len(eq_sys)):
        all_ok = False
        for j in xrange(len(eq_sys)):
            eqj, r_varj, ths_j = eq_sys[j]
            eqi, r_vari, ths_i = eq_sys[i]
            if i == j or set(ths_i) != set(ths_j) or set(r_varj) != set(r_vari):
                continue
            th = ths_i[0]
            C, S = cos(th), sin(th)
            r = r_vari[0]
            X1, Y1, Z1, i_ok = get_coefs(eqi, S, r)
            X2, Y2, Z2, j_ok = get_coefs(eqj, C, r)
            all_ok = j_ok and i_ok and not eqi.has(C) and not eqj.has(S)
            if all_ok:
                eq_type = 5
                break
            X1, Y1, Z1, i_ok = get_coefs(eqi, S, C)
            X2, Y2, Z2, j_ok = get_coefs(eqj, C, S)
            i_ok &= X1.has(r) and not Z1.has(r) and Y1 == ZERO
            j_ok &= X2.has(r) and not Z2.has(r) and Y2 == ZERO
            all_ok = j_ok and i_ok
            if all_ok:
                eq_type = 4
                X1 /= r
                X2 /= r
                break
            else:
                eq_type = -1
        if not all_ok or eq_type == -1:
            continue
        print "#Solving type", eq_type
        if eq_type == 4:
            solve_type_4(symo, X1, Y1, X2, Y2, th, r)
        else:
            solve_type_5(symo, X1, Y1, Z1, X2, Y2, Z2, th, r)
        known |= {th, r}
        return True
    return False

def match_coef(A1,A2,B1,B2):
    return A1 == A2 and B1 == B2 or A1 == -A2 and B1 == -B2

def try_solve_3(symo, eq_sys, known):
    for i in xrange(len(eq_sys)):
        all_ok = False
        for j in xrange(len(eq_sys)):
            eqj, r_varj, ths_i = eq_sys[j]
            eqi, r_vari, ths_j = eq_sys[i]
            if i == j or set(ths_i) != set(ths_j):
                continue
            th1 = ths_i[0]
            th2 = ths_i[1]
            C1, S1 = cos(th1), sin(th1)
            C2, S2 = cos(th2), sin(th2)
            X1, Y1, ZW1, i_ok = get_coefs(eqi, C1, S1)
            X2, Y2, ZW2, j_ok = get_coefs(eqj, S1, C1)
            Y2 = -Y2
            V1, W1, Z1, iw_ok = get_coefs(ZW1, C2, S2)
            V2, W2, Z2, jw_ok = get_coefs(ZW2, S2, C2)
            W2 = -W2
            all_ok = j_ok and i_ok and jw_ok and iw_ok
            if X1 == 0 or Y1 == 0:
                X1, Y1, V1, W1 = V1, W1, X1, Y1
                X2, Y2, V2, W2 = V2, W2, X2, Y2
                th1, th2 = th2, th1
            all_ok &= match_coef(X1,X2,Y1,Y2) and match_coef(V1,V2,W1,W2)
            if W1 == W2 and Y1 == -Y2:
                eps = -1
            else:
                eps = 1
            for a in (X1,X2,Y1,Y2):
                all_ok &= not a.has(C2)
                all_ok &= not a.has(S2)
            if all_ok:
                break
        if not all_ok:
            continue
        print "#Solving type 6, 7"
        solve_type_7(symo, V1, W1, -X1, -Y1, -Z1, -Z2, eps, th1, th2 )
        known |= {th1, th2}
        return True
    return False

def try_solve_4(symo, eq_sys, known):
    res = False
    for i in xrange(len(eq_sys)):
        all_ok = False
        for j in xrange(len(eq_sys)):
            eqj, r_varj, ths_i = eq_sys[j]
            eqi, r_vari, ths_j = eq_sys[i]
            if i == j or set(ths_i) != set(ths_j):
                continue
            th12 = ths_i[0] + ths_i[1]
            if eqi.has(sin(ths_i[0])) or eqi.has(cos(ths_i[0])):
                th1 = ths_i[0]
                th2 = ths_i[1]
            else:
                th1 = ths_i[1]
                th2 = ths_i[0]
            C1, S1 = cos(th1), sin(th1)
            C12, S12 = cos(th12), sin(th12)
            X1, Y1, Z1, i_ok = get_coefs(eqi, C1, C12)
            X2, Y2, Z2, j_ok = get_coefs(eqj, S1, S12)
            all_ok = (X1*Y2 == Y1*X2 and i_ok and j_ok)
            all_ok &= not eqi.has(S1) and not eqi.has(S12)
            all_ok &= not eqj.has(C1) and not eqj.has(C12)
            if all_ok:
                break
        if not all_ok:
            continue
        print "#Solving type 8"
        solve_type_8(symo, X1, Y1, Z1, Z2, th1, th2)
        known |= {th1, th2}
        res = True
    return res

def solve_type_2(symo, X, Y, Z, th):
    """Solution for the equation:
    X*S + Y*C = Z
    """
#    print "#X*sin({0}) + Y*cos({0}) = Z".format(th)
    X = symo.replace(symo.CS12_simp(X), 'X', th)
    Y = symo.replace(symo.CS12_simp(Y), 'Y', th)
    Z = symo.replace(symo.CS12_simp(Z), 'Z', th)
    YPS = var('YPS'+str(th))
    if X == ZERO and Y != ZERO:
        C = symo.replace(Z/Y, 'C', th)
        symo.add_to_dict(YPS, (ONE, - ONE))
        symo.add_to_dict(th, atan2(YPS*sqrt(1-C**2), C))
    elif X != ZERO and Y == ZERO:
        S = symo.replace(Z/X, 'S', th)
        symo.add_to_dict(YPS, (ONE, - ONE))
        symo.add_to_dict(th, atan2(S, YPS*sqrt(1-S**2)))
    elif Z == ZERO:
        symo.add_to_dict(YPS, (ONE, ZERO))
        symo.add_to_dict(th, atan2(-Y, X) + YPS*pi)
    else:
        B = symo.replace(X**2 + Y**2, 'B', th)
        D = symo.replace(B - Z**2, 'D', th)
        symo.add_to_dict(YPS, (ONE, - ONE))
        S = symo.replace((X*Z + YPS * Y * sqrt(D))/B, 'S', th)
        C = symo.replace((Y*Z - YPS * X * sqrt(D))/B, 'C', th)
        symo.add_to_dict(th, atan2(S, C))

def solve_type_3(symo, X1, Y1, Z1, X2, Y2, Z2, th):
    """Solution for the system:
    X1*S + Y1*C = Z1
    X2*S + Y2*C = Z2
    """
#    print "#X1*sin({0}) + Y1*cos({0}) = Z1".format(th)
#    print "#X2*sin({0}) + Y2*cos({0}) = Z2".format(th)
    X1 = symo.replace(symo.CS12_simp(X1), 'X1', th)
    Y1 = symo.replace(symo.CS12_simp(Y1), 'Y1', th)
    Z1 = symo.replace(symo.CS12_simp(Z1), 'Z1', th)
    X2 = symo.replace(symo.CS12_simp(X2), 'X2', th)
    Y2 = symo.replace(symo.CS12_simp(Y2), 'Y2', th)
    Z2 = symo.replace(symo.CS12_simp(Z2), 'Z2', th)
    if X1 == ZERO and Y2 == ZERO:
        symo.add_to_dict(th, atan2(Z2/X2, Z1/Y1))
    elif X2 == ZERO and Y1 == ZERO:
        symo.add_to_dict(th, atan2(Z1/X1, Z2/Y2))
    else:
        D = symo.replace(X1*Y2-X2*Y1, 'D', th)
        C = symo.replace((Z2*X1 - Z1*X2)/D, 'C', th)
        S = symo.replace((Z1*Y2 - Z2*Y1)/D, 'S', th)
        symo.add_to_dict(th, atan2(S, C))

def solve_type_4(symo, X1, Y1, X2, Y2, th, r):
    """Solution for the system:
    X1*S*r = Y1
    X2*C*r = Y2
    """
#    print "X1*sin({0})*{1} = Y1".format(th, r)
#    print "X2*cos({0})*{1} = Y2".format(th, r)
    X1 = symo.replace(symo.CS12_simp(X1), 'X1', th)
    Y1 = symo.replace(symo.CS12_simp(Y1), 'Y1', th)
    X2 = symo.replace(symo.CS12_simp(X2), 'X2', th)
    Y2 = symo.replace(symo.CS12_simp(Y2), 'Y2', th)
    YPS = var('YPS' + r)
    symo.add_to_dict(YPS, (ONE, - ONE))
    symo.add_to_dict(r, YPS*sqrt((Y1/X1)**2 + (Y2/X2)**2))
    symo.add_to_dict(th, atan2(Y1/(X1*r), Y2/(X2*r)))

def solve_type_5(symo, X1, Y1, Z1, X2, Y2, Z2, th, r):
    """Solution for the system:
    X1*S = Y1 + Z1*r
    X2*C = Y2 + Z2*r
    """
#    print "#X1*sin({0}) = Y1 + Z1*{1}".format(th, r)
#    print "#X2*cos({0}) = Y2 + Z2*{1}".format(th, r)
    X1 = symo.replace(symo.CS12_simp(X1), 'X1', th)
    Y1 = symo.replace(symo.CS12_simp(Y1), 'Y1', th)
    Z1 = symo.replace(symo.CS12_simp(Z1), 'Z1', th)
    X2 = symo.replace(symo.CS12_simp(X2), 'X2', th)
    Y2 = symo.replace(symo.CS12_simp(Y2), 'Y2', th)
    Z2 = symo.replace(symo.CS12_simp(Z2), 'Z2', th)
    V1 = symo.replace(Y1/X1, 'V1', r)
    W1 = symo.replace(Z1/X1, 'W1', r)
    V2 = symo.replace(Y2/X2, 'V2', r)
    W2 = symo.replace(Z2/X2, 'W2', r)
    solve_square(W1**2 + W2**2, 2*(V1*W1 + V2*W2), V1**2 + V2**2, r)
    solve_type_3(X1, ZERO, Y1 + Z1*r, ZERO, X2, Y2 + Z2*r)

def solve_type_7(symo, V, W, X, Y, Z1, Z2, eps, th_i, th_j):
    """Solution for the system:
    V1*Cj + W1*Sj = X*Ci + Y*Si + Z1
    eps*(V2*Sj - W2*Cj) = X*Si - Y*Ci + Z2
    """
#    print "#V*cos({0}) + W*sin({0}) = X*cos({1}) + Y*sin({1}) + Z1".format(th_j, th_i)
#    print "#eps*(V*sin({0}) - W*cos({0})) = X*sin({1}) - Y*cos({1}) + Z2".format(th_j, th_i)
    V = symo.replace(symo.CS12_simp(V), 'V', th_i)
    W = symo.replace(symo.CS12_simp(W), 'W', th_i)
    X = symo.replace(symo.CS12_simp(X), 'X', th_i)
    Y = symo.replace(symo.CS12_simp(Y), 'Y', th_i)
    Z1 = symo.replace(symo.CS12_simp(Z1), 'Z1', th_i)
    Z2 = symo.replace(symo.CS12_simp(Z2), 'Z2', th_i)
    B1 = symo.replace(2*(Z1*Y + Z2*X), 'B1', th_i)
    B2 = symo.replace(2*(Z1*X - Z2*Y), 'B2', th_i)
    B3 = symo.replace(V**2 + W**2 - X**2 - Y**2 - Z1**2 - Z2**2, 'B3', th_i)
    solve_type_2(symo, B1, B2, B3, th_i)
    Zi1 = symo.replace(X*cos(th_i) + Y*sin(th_i) + Z1, 'Zi1', th_j)
    Zi2 = symo.replace(X*sin(th_i) - Y*cos(th_i) + Z2, 'Zi2', th_j)
    solve_type_3(symo, W, V, Zi1, eps*V, -eps*W, Zi2, th_j)
#    print_eq(symo, "V1", "X*sin({0}) + Y*cos({0}) + Z1".format(th_i))
#    print_eq(symo, "V2", "X*cos({0}) - Y*sin({0}) + Z2".format(th_i))
#    print_eq(symo, "C", "(V1 - V2)/(2*W2)")
#    print_eq(symo, "S", "(V1 + V2)/(2*W1)")
#    print_eq(symo, th_j, "atan2(S, C)")

def solve_type_8(symo, X, Y, Z1, Z2, th_i, th_j):
    """Solution for the system:
    X*Ci + Y*Cij = Z1
    X*Si + Y*Sij = Z2
    """
#    print "#X*cos({0}) + Y*cos({0} + {1}) = Z1".format(th_i, th_j)
#    print "#X*sin({0}) + Y*sin({0} + {1}) = Z2".format(th_i, th_j)
    X = symo.replace(symo.CS12_simp(X), 'X', th_j)
    Y = symo.replace(symo.CS12_simp(Y), 'Y', th_j)
    Z1 = symo.replace(symo.CS12_simp(Z1), 'Z1', th_j)
    Z2 = symo.replace(symo.CS12_simp(Z2), 'Z2', th_j)
    Cj = symo.replace((Z1**2 + Z2**2 - X**2 - Y**2)/(2*X*Y),'C', th_j)
    YPS = var('YPS' + th_j)
    symo.add_to_dict(YPS, (ONE, - ONE))
    symo.add_to_dict(th_j, atan2(YPS*sqrt(1 - Cj**2), Cj))
    Q1 = symo.replace(X + Y*cos(th_j), 'Q1', th_i)
    Q2 = symo.replace(X + Y*sin(th_j), 'Q2', th_i)
    Den = symo.replace(Q1**2+Q2**2, 'Den', th_i)
    Si = symo.replace((Q1*Z2 - Q2*Z1)/Den, 'S', th_i)
    Ci = symo.replace((Q1*Z1 + Q2*Z2)/Den, 'C', th_i)
    symo.add_to_dict(th_i, atan2(Si,Ci))

def solve_square(symo, A, B, C, x):
    """ solution for the equation:
    A*x**2 + B*x + C = 0
    """
    A = symo.replace(A, 'A', x)
    B = symo.replace(B, 'B', x)
    C = symo.replace(C, 'C', x)
    Delta = symo.repalce(B**2 - 4*A*C, 'Delta', x)
    YPS = var('YPS' + x)
    symo.add_to_dict(YPS, (ONE, - ONE))
    symo.add_to_dict(x, (-B + YPS*sqrt(Delta))/(2*A))

def check_coefs(*args):
    res = True
    for a in args:
        res &= not a.has(sin) and not a.has(cos)
    return res

def match_coefs(X1, X2, Y1, Y2):
    return Y1 == Y2 and X1 == -X2 or Y1 == -Y2 and X1 == X2

def get_coefs(eq, A1, A2):
    eqe = eq.expand()
    X = get_max_coef(eqe, A1)
    eqe = eqe.xreplace({A1:0})
    Y = get_max_coef(eqe, A2)
    Z = eqe.xreplace({A2:0})
    is_ok = not X.has(A2) and not X.has(A1) and not Y.has(A2)
    return X, Y, Z, is_ok


print "Inverse geometric model using Paul method"
T = Matrix([var("s1:4"), var("n1:4"), var("a1:4"), var("p1:4")]).row_join(Matrix([0, 0, 0, 1]))
print T.T
#for i in xrange(3):
#    for j in xrange(3):
#        T[i,j] = EMPTY
f2 = None
def a():
#    robo = Robot.RX90()
    robo = Robot.SR400()
    symo = Symoro()
    loop_solve(Robot.SR400(), symo)
#    igm_Paul(robo, symo, T.T)
    f2 = symo.gen_func('IGM_generated',robo.get_q_vec(),robo.get_q_active())

import profile
from timeit import timeit
#print timeit(a,number = 1)
a()
#print timeit(f2,number = 1)
#profile.run('a()', sort = 'cumtime')

