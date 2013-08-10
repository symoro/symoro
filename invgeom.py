"""
This module of SYMORO package provides symbolic
solutions for inverse geompetric problem.

The core symbolic library is sympy.
Needed modules : symoro.py, geometry.py

ECN - ARIA1 2013
"""

#get all transformations
#for each equation try to mach
#with the pattern equation.
#if something is matched, output the solution,
#put the found variables to the constant list

from sympy import Matrix, var, Symbol, sin, cos
from symoro import Symoro, Robot, ZERO
from symoro import get_trig_couple_names, get_max_coef
from geometry import dgm

EMPTY = var("EMPTY")

#dictionary for equation type classification
eq_dict = {(1, 0, 0):0, (0, 1, 0):1, (1, 1, 0):2,
                  (0, 2, 0):3, (0, 2, 1):4}
name_dict = ("type 1", "type 2, 3", "type 4, 5", "type 6, 7", "type 8")

def print_eq(A, B):
    print "    " + str(A) + " = " + str(B)

def igm_Paul(robo, symo, T_ref):
    #get the terminal frames
    terminal = list(set(range(robo.NL))-set(robo.ant))
    for n in terminal:   #terminal:
        print "#solution for the frame", robo.num[n]
        chain = robo.chain(n)
        r_all = set(robo.r[i] for i in chain if robo.sigma[i] == 1)
        th_all = set(robo.theta[i] for i in chain if robo.sigma[i] == 0)
        known = set()
        iTn = dgm(robo, symo, -1, n, fast_form = False,
                  all_trans = True, trig_subs = False)
        iT0 = dgm(robo, symo, n, -1, fast_form = False,
                  all_trans = True, trig_subs = False)
        while True:
            repeat_search = False
            for i in reversed(chain + [-1]):
#                print iT0[(i, -1)]*T_ref
#                print iTn[(i, n)]
#                print "#transformations", i
                M_eq = iT0[(i, -1)]*T_ref - iTn[(i, n)]
#                print M_eq
                while True:
                    cont_search = False
                    eq_candidates = [list() for list_index in range(5)]
                    for e1 in range(3):
                        for e2 in range(4):
                            eq = M_eq[e1, e2]#.expand()
                            if eq.has(EMPTY):
                                continue
#                            print i
#                            print eq
                            names_c, names_s = get_trig_couple_names(eq)
                            th_vars = (eq.atoms(Symbol) & th_all) - known
                            if th_vars:
                                arg_sum = max(at.count_ops()-1 for at in eq.atoms(sin, cos)
                                            if not (at.atoms(Symbol) & th_all) <= known)
                            else:
                                arg_sum = 0
                            r_vars = (eq.atoms(Symbol) & r_all) - known
    #                        if len(r_vars) + len(th_names) < 3:
    #                            print eq
                            eq_features = (len(r_vars), len(th_vars), arg_sum)
                            if eq_dict.has_key(eq_features):
                                eq_key = eq_dict[eq_features]
                                eq_pack = (eq, list(r_vars), list(th_vars))
                                eq_candidates[eq_key].append(eq_pack)
#                    for key_print in range(5):
#                        print key_print
#                        for equation in eq_candidates[key_print]:
#                            print equation
                    cont_search |= try_solve_0(eq_candidates[0], known)
                    cont_search |= try_solve_1(eq_candidates[1], known)
                    cont_search |= try_solve_2(eq_candidates[2] + eq_candidates[1], known)
                    cont_search |= try_solve_3(eq_candidates[3], known)
                    cont_search |= try_solve_4(eq_candidates[4], known)
                    repeat_search |= cont_search
                    if not cont_search or th_all | r_all <= known:
                        break
                if th_all | r_all <= known:
                    break
            if not repeat_search or th_all | r_all <= known:
                break

def try_solve_0(eq_sys, known):
    res = False
    for eq, [r], th_names in eq_sys:
        X = get_max_coef(eq, r)
        if X != 0:
            Y = X*r - eq
            print "type 1"
            print_eq("X", X)
            print_eq("Y", Y)
            print_eq(r, "Y/X")
            known.add(r)
            res = True
    return res

def try_solve_1(eq_sys, known):
    res = False
    for i in range(len(eq_sys)):
        eqi, r_vari, [th_i ] = eq_sys[i]
        if th_i  in known:
            continue
        Xi, Yi, Zi, i_ok = get_coefs(eqi, sin(th_i ), cos(th_i ))
        i_ok &= sum([Xi == ZERO, Yi == ZERO, Zi == ZERO]) <= 1
        if not i_ok:
            continue
        j_ok = False
        for j in range(i+1, len(eq_sys)):
            eqj, r_varj, [th_j] = eq_sys[j]
            if th_i  == th_j:
                Xj, Yj, Zj, j_ok = get_coefs(eqj, sin(th_j), cos(th_j))
                j_ok &= (Xi*Yj != Xj*Yi)
                if j_ok:
                    break
        if j_ok:
            print "#Solving type 3"
            solve_type_3(Xi, Yi, -Zi, Xj, Yj, -Zj, th_i)
        else:
            print "#Solving type 2"
            solve_type_2(Xi, Yi, -Zi, th_i)
        known.add(th_i )
        res = True
    return res

def try_solve_2(eq_sys, known):
    if all(len(r_var) == 0 for eq, r_var, ths in eq_sys):
        return False
    for i in range(len(eq_sys)):
        all_ok = False
        for j in range(len(eq_sys)):
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
            solve_type_4(X1, Y1, X2, Y2, th, r)
        else:
            solve_type_5(X1, Y1, Z1, X2, Y2, Z2, th, r)
        known |= {th, r}
        return True
    return False

def try_solve_3(eq_sys, known):
    for i in range(len(eq_sys)):
        all_ok = False
        for j in range(len(eq_sys)):
            eqj, r_varj, ths_i = eq_sys[j]
            eqi, r_vari, ths_j = eq_sys[i]
            if i == j or set(ths_i) != set(ths_j):
                continue
            th1 = ths_i[0]
            th2 = ths_i[1]
            C1, S1 = cos(th1), sin(th1)
            C2, S2 = cos(th2), sin(th2)
            def match_coef(A1,A2,B1,B2):
                return A1 == A2 and B1 == B2 or A1 == -A2 and B1 == -B2
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
        print "#",eqi
        print "#",eqj
        solve_type_7(V1, W1, -X1, -Y1, -Z1, -Z2, eps, th1, th2 )
        known |= {th1, th2}
        return True
    return False

def try_solve_4(eq_sys, known):
    res = False
    for i in range(len(eq_sys)):
        all_ok = False
        for j in range(len(eq_sys)):
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
        solve_type_8(X1, Y1, Z1, Z2, th1, th2)
        known |= {th1, th2}
        res = True
    return res

def solve_type_2(X, Y, Z, th):
    """Solution for the equation:
    X*S + Y*C = Z
    """
    print "#X*sin({0}) + Y*cos({0}) = Z".format(th)
    print_eq("X", X)
    print_eq("Y", Y)
    print_eq("Z", Z)
    if X == ZERO and Y != ZERO:
        print_eq("C".format(th), "Z/Y")
        print_eq(th, "atan2(sqrt(1-C**2), C)")
        print "or"
        print_eq(th, "atan2(-sqrt(1-C**2), C)")
    elif X != ZERO and Y == ZERO:
        print_eq("S".format(th), "Z/X")
        print_eq(th, "atan2(S, sqrt(1-S**2))")
        print "or"
        print_eq(th, "atan2(S, -sqrt(1-S**2))")
    elif Z == ZERO:
        print_eq(th, "atan2(-Y, X)")
        print "or"
        print_eq(th, "atan2(-Y, X) + pi")
    else:
        print_eq("B", "X**2 + Y**2")
        print_eq("D", "B - Z**2")
        print_eq("S", "(X*Z + Y * sqrt(D))/B")
        print_eq("C", "(Y*Z - X * sqrt(D))/B")
        print "or"
        print_eq("S", "(X*Z - Y * sqrt(D))/B")
        print_eq("C", "(Y*Z + X * sqrt(D))/B")
        print_eq(th, "atan2(S, C)")

def solve_type_3(X1, Y1, Z1, X2, Y2, Z2, th):
    """Solution for the system:
    X1*S + Y1*C = Z1
    X2*S + Y2*C = Z2
    """
    print "#X1*sin({0}) + Y1*cos({0}) = Z1".format(th)
    print "#X2*sin({0}) + Y2*cos({0}) = Z2".format(th)
    print_eq("X1", X1)
    print_eq("Y1", Y1)
    print_eq("Z1", Z1)
    print_eq("X2", X2)
    print_eq("Y2", Y2)
    print_eq("Z2", Z2)
    if X1 == ZERO and Y2 == ZERO:
        print_eq(th, "atan2(Z2/X2, Z1/Y1)")
    elif X2 == ZERO and Y1 == ZERO:
        print_eq(th, "atan2(Z1/X1, Z2/Y2)")
    else:
        print_eq("D", "(X1*Y2-X2*Y1)")
        print_eq("C", "(Z2*X1 - Z1*X2)/D")
        print_eq("S", "(Z1*Y2 - Z2*Y1)/D")
        print_eq(th, "atan2(S, C)")

def solve_type_4(X1, Y1, X2, Y2, th, r):
    """Solution for the system:
    X1*S*r = Y1
    X2*C*r = Y2
    """
    print "X1*sin({0})*{1} = Y1".format(th, r)
    print "X2*cos({0})*{1} = Y2".format(th, r)
    print_eq("X1", X1)
    print_eq("Y1", Y1)
    print_eq("X2", X2)
    print_eq("Y2", Y2)
    print_eq(r, "sqrt((Y1/X1)**2 + (Y2/X2)**2)")
    print "or"
    print_eq(r, "-sqrt((Y1/X1)**2 + (Y2/X2)**2)")
    print_eq(th, "atan2((Y1/(X1*{0}), Y2/(X2*{0}))".format(r))

def solve_type_5(X1, Y1, Z1, X2, Y2, Z2, th, r):
    """Solution for the system:
    X1*S = Y1 + Z1*r
    X2*C = Y2 + Z2*r
    """
    print "#X1*sin({0}) = Y1 + Z1*{1}".format(th, r)
    print "#X2*cos({0}) = Y2 + Z2*{1}".format(th, r)
    print_eq("X1", X1)
    print_eq("Y1", Y1)
    print_eq("Z1", Z1)
    print_eq("X2", X2)
    print_eq("Y2", Y2)
    print_eq("Z2", Z2)
    print_eq("V1", "Y1/X1")
    print_eq("W1", "Z1/X1")
    print_eq("V2", "Y2/X2")
    print_eq("W2", "Z2/X2")
    solve_square("W1**2 + W2**2", "2*(V1*W1 + V2*W2)", "V1**2 + V2**2", r)
    solve_type_3(X1, ZERO, Y1 + Z1*r, ZERO, X2, Y2 + Z2*r)

def solve_type_7(V, W, X, Y, Z1, Z2, eps, th_i, th_j):
    """Solution for the system:
    V1*Cj + W1*Sj = X*Ci + Y*Si + Z1
    eps*(V2*Sj - W2*Cj) = X*Si - Y*Ci + Z2
    """
    print "#V*cos({0}) + W*sin({0}) = X*cos({1}) + Y*sin({1}) + Z1".format(th_j, th_i)
    print "#eps*(V*sin({0}) - W*cos({0})) = X*sin({1}) - Y*cos({1}) + Z2".format(th_j, th_i)
    print_eq("#V", V)
    print_eq("#W", W)
    print_eq("#eps",eps)
    print_eq("#X", X)
    print_eq("#Y", Y)
    print_eq("#Z1", Z1)
    print_eq("#Z2", Z2)
    B1 = 2*(Z1*Y + Z2*X)
    B2 = 2*(Z1*X - Z2*Y)
    B3 = V**2 + W**2 - X**2 - Y**2 - Z1**2 - Z2**2
    solve_type_2(B1, B2, B3, th_i)
    Zi1 = X*cos(th_i) + Y*sin(th_i) + Z1
    Zi2 = X*sin(th_i) - Y*cos(th_i) + Z2
    solve_type_3(W, V, Zi1, eps*V, -eps*W, Zi2, th_j)
#    print_eq("V1", "X*sin({0}) + Y*cos({0}) + Z1".format(th_i))
#    print_eq("V2", "X*cos({0}) - Y*sin({0}) + Z2".format(th_i))
#    print_eq("C", "(V1 - V2)/(2*W2)")
#    print_eq("S", "(V1 + V2)/(2*W1)")
#    print_eq(th_j, "atan2(S, C)")

def solve_type_8(X, Y, Z1, Z2, thi, thj):
    """Solution for the system:
    X*Ci + Y*Cij = Z1
    X*Si + Y*Sij = Z2
    """
    print "#X*cos({0}) + Y*cos({0} + {1}) = Z1".format(thi, thj)
    print "#X*sin({0}) + Y*sin({0} + {1}) = Z2".format(thi, thj)
    print_eq("X", X)
    print_eq("Y", Y)
    print_eq("Z1", Z1)
    print_eq("Z2", Z2)
    print_eq("Cj","(Z1**2 + Z2**2 - X**2 - Y**2)/(2*X*Y)")
    print_eq(thj, "atan2(sqrt(1 - Cj**2), Cj)")
    print "or"
    print_eq(thj, "atan2(-sqrt(1 - Cj**2), Cj)")
    print_eq("B1", "X + Y*cos({0})".format(thj))
    print_eq("B2", "Y*sin({0})".format(thj))
    print_eq("Si", "(B1*Z2 - B2*Z1)/(B1**2+B2**2)")
    print_eq("Ci", "(B1*Z1 + B2*Z2)/(B1**2+B2**2)")
    print_eq(thi, "atan2(Si,Ci)")

def solve_square(A, B, C, x):
    """ solution for the equation:
    A*x**2 + B*x + C = 0
    """
    print_eq("A", A)
    print_eq("B", B)
    print_eq("C", C)
    print_eq("Delta", "B**2 - 4*A*C")
    print_eq(x, "(-B + sqrt(Delta))/(2*A)")
    print "or"
    print_eq(x, "(-B - sqrt(Delta))/(2*A)")

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
    eqe = eqe.subs({A1:0})
    Y = get_max_coef(eqe, A2)
    Z = eqe.subs({A2:0})
    is_ok = not X.has(A2) and not X.has(A1) and not Y.has(A2)
    return X, Y, Z, is_ok

print "Inverse geometric model using Paul method"
T = Matrix([var("s1:4"), var("n1:4"), var("a1:4"), var("p1:4")]).row_join(Matrix([0, 0, 0, 1]))
print T.T
#for i in range(3):
#    for j in range(3):
#        T[i,j] = EMPTY
igm_Paul(Robot.RX90(), Symoro(), T.T)

