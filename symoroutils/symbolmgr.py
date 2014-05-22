# -*- coding: utf-8 -*-


"""This module contains the Symbol Manager tools."""

import itertools

from sympy import sin, cos
from sympy import Symbol, Matrix, Expr
from sympy import Mul, Add, factor, var, sympify

from symoroutils import filemgr
from symoroutils import tools


class SymbolManager(object):
    """Symbol manager, responsible for symbol replacing, file writing."""
    def __init__(self, file_out='disp', sydi=dict()):
        """Default values correspond to empty dictionary and screen output.
        """
        self.file_out = file_out
        """Output descriptor. Can be None, 'disp', file
        defines the output destination"""
        self.sydi = dict((k, sydi[k]) for k in sydi)
        """Dictionary. All the substitutions are saved in it"""
        self.revdi = dict((sydi[k], k) for k in sydi)
        """Dictionary. Revers to the self.sydi"""
        self.order_list = sydi.keys()
        """keeps the order of variables to be compute"""

    def simp(self, sym):
        sym = factor(sym)
        new_sym = tools.ONE
        for expr in Mul.make_args(sym):
            if expr.is_Pow:
                expr, pow_val = expr.args
            else:
                pow_val = 1
            expr = self.C2S2_simp(expr)
            expr = self.CS12_simp(expr, silent=True)
            new_sym *= expr**pow_val
        return new_sym

    def C2S2_simp(self, sym):
        """
        Example
        =======
        >> print C2S2_simp(sympify("-C**2*RL + S*(D - RL*S)"))
        D*S - RL
        """
        if not sym.is_Add:
            repl_dict = {}
            for term in sym.atoms(Add):
                repl_dict[term] = self.C2S2_simp(term)
            sym = sym.xreplace(repl_dict)
            return sym
        names, short_form = tools.trignometric_info(sym)
        for name in names:
            if short_form:
                cos_term, sin_term = tools.cos_sin_syms(name)
            else:
                cos_term, sin_term = cos(name), sin(name)
            sym = self.try_opt(
                tools.ONE, None, sin_term**2, cos_term**2, sym
            )
        return sym

    def CS12_simp(self, sym, silent=False):
        """
        Example
        =======
        >> print SymbolManager().CS12_simp(sympify("C2*C3 - S2*S3"))
        C23 = C2*C3 - S2*S3
        C23
        >> print SymbolManager().CS12_simp(sympify("C2*S3*R + S2*C3*R"))
        S23 = C2*S3 + S2*C3
        R*S23
        """
        if not sym.is_Add:
            repl_dict = {}
            for term in sym.atoms(Add):
                repl_dict[term] = self.CS12_simp(term)
            sym = sym.xreplace(repl_dict)
            return sym
        names, short_form = tools.trignometric_info(sym)
        names = list(names)
        if short_form:
            names.sort()
        sym2 = sym
        for n1, n2 in itertools.combinations(names, 2):
            if short_form:
                C1, S1 = tools.cos_sin_syms(n1)
                C2, S2 = tools.cos_sin_syms(n2)
                np1, nm1 = tools.get_pos_neg(n1)
                np2, nm2 = tools.get_pos_neg(n2)
                n12 = tools.ang_sum(np1, np2, nm1, nm2)
                nm12 = tools.ang_sum(np1, nm2, nm1, np2)
                C12, S12 = tools.cos_sin_syms(n12)
                C1m2, S1m2 = tools.cos_sin_syms(nm12)
            else:
                C1, S1 = cos(n1), sin(n1)
                C2, S2 = cos(n2), sin(n2)
                C12, S12 = cos(n1+n2), sin(n1+n2)
                C1m2, S1m2 = cos(n1-n2), sin(n1-n2)
            sym2 = self.try_opt(S12, S1m2, S1*C2, C1*S2, sym2, silent)
            sym2 = self.try_opt(C12, C1m2, C1*C2, -S1*S2, sym2, silent)
        if sym2 != sym:
            return self.CS12_simp(sym2, silent)
        else:
            return sym

    def try_opt(self, A, Am, B, C, old_sym, silent=False):
        """Replaces B + C by A or B - C by Am.
        Chooses the best option.
        """
        Bcfs = tools.get_max_coef_list(old_sym, B)
        Ccfs = tools.get_max_coef_list(old_sym, C)
        if Bcfs != [] and Ccfs != []:
            Res = old_sym
            Res_tmp = Res
            for coef in Bcfs:
                Res_tmp += A*coef - B*coef - C*coef
                if tools.sym_less(Res_tmp, Res):
                    Res = Res_tmp
            if tools.sym_less(Res, old_sym) and Am is None:
                if not A.is_number and not silent:
                    self.add_to_dict(A, B + C)
                return Res
            elif Am is not None:
                Res2 = old_sym
                Res_tmp = Res2
                for coef in Bcfs:
                    Res_tmp += Am*coef - B*coef + C*coef
                    if tools.sym_less(Res_tmp, Res2):
                        Res2 = Res_tmp
                if tools.sym_less(Res2, Res) and tools.sym_less(Res2, old_sym):
                    if not Am.is_number and not silent:
                        self.add_to_dict(Am, B - C)
                    return Res2
                elif tools.sym_less(Res, old_sym):
                    if not A.is_number and not silent:
                        self.add_to_dict(A, B + C)
                    return Res
        return old_sym

    def add_to_dict(self, new_sym, old_sym):
        """Internal function.
        Extends symbol dictionary by (new_sym, old_sym) pair
        """
        new_sym = sympify(new_sym)
        if new_sym.as_coeff_Mul()[0] == -tools.ONE:
            new_sym = -new_sym
            old_sym = -old_sym
        if new_sym not in self.sydi:
            self.sydi[new_sym] = old_sym
            self.revdi[old_sym] = new_sym
            self.order_list.append(new_sym)
            self.write_equation(new_sym, old_sym)

    def trig_replace(self, M, angle, name):
        """Replaces trigonometric expressions cos(x)
        and sin(x) by CX and SX

        Parameters
        ==========
        M: var or Matrix
            Object of substitution
        angle: var
            symbol that stands for the angle value
        name: int or string
            brief name X for the angle

        Notes
        =====
        The cos(x) and sin(x) will be replaced by CX and SX,
        where X is the name and x is the angle
        """
        if not isinstance(angle, Expr) or angle.is_number:
            return M
        cos_sym, sin_sym = tools.cos_sin_syms(name)
        sym_list = [(cos_sym, cos(angle)), (sin_sym, sin(angle))]
        subs_dict = {}
        for sym, sym_old in sym_list:
            if sym_old.has(-1):
                subs_dict[-sym_old] = -sym
            else:
                subs_dict[sym_old] = sym
            self.add_to_dict(sym, sym_old)
        for i1 in xrange(M.shape[0]):
            for i2 in xrange(M.shape[1]):
                M[i1, i2] = M[i1, i2].subs(subs_dict)
        return M

    def replace(self, old_sym, name, index='', forced=False):
        """Creates a new symbol for the symbolic expression old_sym.

        Parameters
        ==========
        old_sym: var
            Symbolic expression to be substituted
        name: string or var
            denotion of the expression
        index: int or string, optional
            will be attached to the name. Usualy used for link or joint number.
            Parameter exists for usage convenience
        forced: bool, optional
            If True, the new symbol will be created even if old symbol
            is a simple expression

        Notes
        =====
        Generaly only complex expressions, which contain + - * / ** operations
        will be replaced by a new symbol
        """
        if not forced:
            if not isinstance(old_sym, Expr):
                return old_sym
            inv_sym = -old_sym
            if old_sym.is_Atom or inv_sym.is_Atom:
                return old_sym
            for i in (1, -1):
                if i * old_sym in self.revdi:
                    return i * self.revdi[i * old_sym]
        new_sym = var(str(name) + str(index))
        self.add_to_dict(new_sym, old_sym)
        return new_sym

    def mat_replace(self, M, name, index='',
                    forced=False, skip=0, symmet=False):
        """Replaces each element in M by symbol

        Parameters
        ==========
        M: Matrix
            Object of substitution
        name: string
            denotion of the expression
        index: int or string, optional
            will be attached to the name. Usualy used for link
            or joint number. Parameter exists for usage convenience
        forced: bool, optional
            If True, the new symbol will be created even if old symbol
            is a simple expression
        skip: int, optional
            Number of bottom rows of the matrix, which will be skipped.
            Used in case of Transformation matrix and forced = True.
        symmet: bool, optional
            If true, only for upper triangle part of the matrix
            symbols will be created. The bottom triangle part the
            same symbols will be used


        Returns
        =======
        M: Matrix
            Matrix with all the elements replaced

        Notes
        =====
        -Each element M_ij will be replaced by
            symbol name + i + j + index
        -There are two ways to use this function (examples):
            1)  >>> A = B+C+...
                >>> symo.mat_replace(A, 'A')
                # for the case when expression B+C+... is too big
            2)  >>> A = symo.mat_replace(B+C+..., 'A')
                # for the case when B+C+... is small enough
        """
        for i2 in xrange(M.shape[1]):
            for i1 in xrange(M.shape[0] - skip):
                if symmet and i2 < i1:
                    M[i1, i2] = M[i2, i1]
                    continue
                if M.shape[1] > 1:
                    name_index = name + str(i1 + 1) + str(i2 + 1)
                else:
                    name_index = name + str(i1 + 1)
                M[i1, i2] = self.replace(M[i1, i2], name_index, index, forced)
        return M

    def unfold(self, expr):
        """Unfold the expression using the dictionary.

        Parameters
        ==========
        expr: symbolic expression
            Symbolic expression to be unfolded

        Returns
        =======
        expr: symbolic expression
            Unfolded expression
        """
        while self.sydi.keys() & expr.atoms():
            expr = expr.subs(self.sydi)
        return expr

    def write_param(self, name, header, robo, N):
        """Low-level function for writing the parameters table

        Parameters
        ==========
        name: string
            the name of the table
        header: list
            the table header
        robo: Robot
            Instance of parameter container
        N: list of int
            Indices for which parameter rows will be written
        """
        self.write_line(name)
        self.write_line(tools.l2str(header))
        for j in N:
            params = robo.get_param_vec(header, j)
            self.write_line(tools.l2str(params))
        self.write_line()

    def write_params_table(self, robo, title='', geom=True, inert=False,
                           dynam=False, equations=True,
                           inert_name='Dynamic inertia parameters'):
        """Writes the geometric parameters table

        Parameters
        ==========
        robo: Robot
            Instance of the parameter container.
        title: string
            The document title.

        Notes
        =====
        The synamic model generation program can be started with this function
        """
        if title != '':
            self.write_line(title)
            self.write_line()
        if geom:
            self.write_param('Geometric parameters', robo.get_geom_head(),
                             robo, range(1, robo.NF))
        if inert:
            if robo.is_mobile:
                start_frame = 0
            else:
                start_frame = 1

            self.write_param(inert_name, robo.get_dynam_head(),
                             robo, range(start_frame, robo.NL))
        if dynam:
            self.write_param('External forces and joint parameters',
                             robo.get_ext_dynam_head(),
                             robo, range(1, robo.NL))
            self.write_param('Base velicities parameters',
                             robo.get_base_vel_head(),
                             robo, [0, 1, 2])
        if equations:
            self.write_line('Equations:')

    def unknown_sep(self, eq, known):
        """If there is a sum inside trigonometric function and
        the atoms are not the subset of 'known',
        this function will replace the trigonometric symbol bu sum,
        trying to separate known and unknown terms
        """
        if not isinstance(eq, Expr) or eq.is_number:
            return eq
        while True:
            res = False
            trigs = eq.atoms(sin, cos)
            for trig in trigs:
                args = trig.args[0].atoms()
                if args & known and not args <= known and trig in self.sydi:
                    eq = eq.subs(trig, self.sydi[trig]).expand()
                    res = True
            if not res:
                break
        return eq

    def write_equation(self, A, B):
        """Writes the equation A = B into the output

        Parameters
        ==========
        A: expression or var
            left-hand side of the equation.
        B: expression or var
            right-hand side of the equation
        """
        self.write_line(str(A) + ' = ' + str(B))

    def write_line(self, line=''):
        """Writes string data into tha output with new line symbol

        Parameters
        ==========
        line: string, optional
            Data to be written. If empty, it adds an empty line
        """
        if self.file_out == 'disp':
            print line
        elif self.file_out is not None:
            self.file_out.write(str(line) + '\n')

    def file_open(self, robo, ext):
        """
        Initialize file stream

        Parameters
        ==========
        robo: Robot instance
            provides the robot's name
        ext: string
            provides the file name extention
        """
        fname = filemgr.make_file_path(robo, ext)
        self.file_out = open(fname, 'w')

    def file_close(self):
        """
        Initialize file stream

        Parameters
        ==========
        robo: Robot instance
            provides the robot's name
        ext: string
            provides the file name extention
        """
        if self.file_out is not None:
            self.write_line('*=*')
            self.file_out.close()

    def gen_fheader(self, name, *args):
        fun_head = []
        fun_head.append('def %s_func(*args):\n' % name)
        imp_s_1 = 'from numpy import pi, sin, cos, sign\n'
        imp_s_2 = 'from numpy import array, arctan2 as atan2, sqrt\n'
        fun_head.append('    %s' % imp_s_1)
        fun_head.append('    %s' % imp_s_2)
        for i, var_list in enumerate(args):
            v_str_list = self.convert_syms(args[i], True)
            fun_head.append('    %s=args[%s]\n' % (v_str_list, i))
        return fun_head

    def convert_syms(self, syms, rpl_liter=False):
        """Converts 'syms' structure to sintactically correct string

        Parameters
        ==========
        syms: list, Matrix or tuple of them
        rpl_liter: bool
            if true, all literals will be replaced with _
            It is done to evoid expression like [x, 0] = args[1]
            Because it will cause exception of assigning to literal
        """
        if isinstance(syms, tuple) or isinstance(syms, list):
            syms = [self.convert_syms(item, rpl_liter) for item in syms]
            res = '['
            for i, s in enumerate(syms):
                res += s
                if i < len(syms) - 1:
                    res += ','
            res += ']'
            return res
        elif isinstance(syms, Matrix):
            res = '['
            for i in xrange(syms.shape[0]):
                res += self.convert_syms(list(syms[i, :]), rpl_liter)
                if i < syms.shape[0] - 1:
                    res += ','
            res += ']'
            return res
        else:
            if rpl_liter and sympify(syms).is_number:
                return '_'
            else:
                return str(syms)

    def extract_syms(self, syms):
        """ returns set of all symbols from list or matrix
        or tuple of them
        """
        if isinstance(syms, tuple) or isinstance(syms, list):
            atoms = (self.extract_syms(item) for item in syms)
            return reduce(set.__or__, atoms, set())
        elif isinstance(syms, Matrix):
            return self.extract_syms(list(syms))
        elif isinstance(syms, Expr):
            return syms.atoms(Symbol)
        else:
            return set()

    def sift_syms(self, rq_syms, wr_syms):
        """Returns ordered list of variables to be compute
        """
        order_list = []   # vars that are defined in sydi
        for s in reversed(self.order_list):
            if s in rq_syms and not s in wr_syms:
                order_list.insert(0, s)
                s_val = self.sydi[s]
                if isinstance(s_val, Expr):
                    atoms = s_val.atoms(Symbol)
                    rq_syms |= {s for s in atoms if not s.is_number}
        rq_vals = [s for s in rq_syms if not (s in self.sydi or s in wr_syms)]
            # required vars that are not defined in sydi
            # will be set to '1.'
        return rq_vals + order_list

    def gen_fbody(self, name, to_return, wr_syms, multival):
        """Generates list of string statements of the function that
        computes symbolf from to_return.  wr_syms are considered to
        be known
        """
        # final symbols to be compute
        syms = self.extract_syms(to_return)
        # defines order of computation
        order_list = self.sift_syms(syms, wr_syms)
        # list of instructions in final function
        fun_body = []
        # will be switched to true when branching detected
        space = '    '
        folded = 1    # indentation = 1 + number of 'for' statements

        for s in order_list:
            if s not in self.sydi:
                item = '%s%s=1.\n' % (space * folded, s)
            elif isinstance(self.sydi[s], tuple):
                multival = True
                item = '%sfor %s in %s:\n' % (space * folded, s, self.sydi[s])
                folded += 1
            else:
                item = '%s%s=%s\n' % (space * folded, s, self.sydi[s])
            fun_body.append(item)
        ret_expr = self.convert_syms(to_return)
        if multival:
            fun_body.insert(0, '    %s_result=[]\n' % (name))
            item = '%s%s_result.append(%s)\n' % (space*folded, name, ret_expr)
        else:
            item = '    %s_result=%s\n' % (name, ret_expr)
        fun_body.append(item)
        fun_body.append('    return %s_result\n' % (name))
        return fun_body

    def gen_func(self, name, to_return, args, multival=False):
        """ Returns function that computes what is in to_return
        using *args as arguments

         Parameters
        ==========
        name: string
            Future function's name, must be different for
            different fucntions
        to_return: list, Matrix or tuple of them
            Determins the shape of the output and symbols inside it
        *args: any number of lists, Matrices or tuples of them
            Determins the shape of the input and symbols
            names to assigned

        Notes
        =====
        -All unassigned used symbols will be set to '1.0'.
        -This function must be called only after the model that
            computes symbols in to_return have been generated.
        """
        fun_head = self.gen_fheader(name, args)
        wr_syms = self.extract_syms(args)   # set of defined symbols
        fun_body = self.gen_fbody(name, to_return, wr_syms, multival)
        fun_string = "".join(fun_head + fun_body)
        exec fun_string
        return eval('%s_func' % name)


