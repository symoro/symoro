# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""This module contains the Symbol Manager tools."""


import os

from sympy import sin, cos
from sympy import Symbol, Matrix, Expr
from sympy import var, sympify

from symoroutils import filemgr
from symoroutils import tools
from genfunc import gen_fheader_matlab, gen_fbody_matlab

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

    def clear(self):
        self.revdi.clear()
        self.sydi.clear()
        self.order_list = []

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

    def trig_replace(self, M):
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
        all_cos = M.atoms(cos)
        all_sin = M.atoms(sin)
        subs_dict = {}

        def create_var(f):
            if isinstance(f, cos):
                name = 'C'
            elif isinstance(f, sin):
                name = 'S'
            syms = f.atoms(Symbol)
            for s in syms:
                name += '_%s' % s
            return var(name)

        for func in all_cos | all_sin:
            if func not in self.revdi:
                subs_dict[func] = create_var(func)
                self.add_to_dict(subs_dict[func], func)
            else:
                subs_dict[func] = self.revdi[func]
        for i1 in xrange(M.shape[0]):
            for i2 in xrange(M.shape[1]):
                M[i1, i2] = M[i1, i2].subs(subs_dict)
        return M

    #TODO remove index
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
        if M.shape[0] > 9:
            form2 = '%02d%02d'
        else:
            form2 = '%d%d'
        for i2 in xrange(M.shape[1]):
            for i1 in xrange(M.shape[0] - skip):
                if symmet and i1 < i2:
                    M[i1, i2] = M[i2, i1]
                    continue
                if M.shape[1] > 1:
                    name_index = name + form2 % (i1 + 1, i2 + 1)
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
        while set(self.sydi.keys()) & expr.atoms():
            expr = expr.subs(self.sydi)
        return expr

    def mat_unfold(self, mat):
        for i in xrange(mat.shape[0]):
            for j in xrange(mat.shape[1]):
                if isinstance(mat[i, j], Expr):
                    mat[i, j] = self.unfold(mat[i, j])
        return mat

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
            if robo.is_floating or robo.is_mobile:
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

#    def unknown_sep(self, eq, known):
#        """If there is a sum inside trigonometric function and
#        the atoms are not the subset of 'known',
#        this function will replace the trigonometric symbol bu sum,
#        trying to separate known and unknown terms
#        """
#        if not isinstance(eq, Expr) or eq.is_number:
#            return eq
#        while True:
#            res = False
#            trigs = eq.atoms(sin, cos)
#            for trig in trigs:
#                args = trig.args[0].atoms()
#                if args & known and not args <= known and trig in self.sydi:
#                    eq = eq.subs(trig, self.sydi[trig]).expand()
#                    res = True
#            if not res:
#                break
#        return eq

    def write_equation(self, A, B):
        """Writes the equation A = B into the output

        Parameters
        ==========
        A: expression or var
            left-hand side of the equation.
        B: expression or var
            right-hand side of the equation
        """
        self.write_line(str(A) + ' = ' + str(B) + ';')

    def write_line(self, line=''):
        """Writes string data into tha output with new line symbol

        Parameters
        ==========
        line: string, optional
            Data to be written. If empty, it adds an empty line
        """
        if self.file_out == 'disp':
            print(line)
        elif self.file_out is not None:
            self.file_out.write(str(line) + '\n')

    def flushout(self):
        """
        Flush the buffer and make sure the data is written to the disk
        """
        self.file_out.flush()
        if self.file_out != 'disp':
            os.fsync(self.file_out.fileno())

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
        fname = filemgr.get_file_path(robo, ext)
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
        fun_head.append('def %s(*args):\n' % name)
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
        elif rpl_liter and sympify(syms).is_number:
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

    def gen_fbody(self, name, to_return, args):
        """Generates list of string statements of the function that
        computes symbolf from to_return.  wr_syms are considered to
        be known
        """
         # set of defined symbols
        wr_syms = self.extract_syms(args)
        # final symbols to be compute
        syms = self.extract_syms(to_return)
        # defines order of computation
        order_list = self.sift_syms(syms, wr_syms)
        # list of instructions in final function
        fun_body = []
        # will be switched to true when branching detected
        space = '    '
        folded = 1    # indentation = 1 + number of 'for' statements
        multival = False
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

    def gen_func_string(self, name, to_return, args, syntax='python'):
        #TODO self, name, toret, *args, **kwargs
        """ Returns function string. The rest is the same as for
        gen_func

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
        #if kwargs.get
        if syntax == 'python':
            fun_head = self.gen_fheader(name, args)
            fun_body = self.gen_fbody(name, to_return, args)
        elif syntax == 'matlab':
            fun_head = gen_fheader_matlab(self, name, args, to_return)
            fun_body = gen_fbody_matlab(self, name, to_return, args)
        fun_string = "".join(fun_head + fun_body)
        return fun_string

    def gen_func(self, name, to_return, args):
        """ Returns function that computes what is in to_return
        using args as arguments

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
        s = self.gen_func_string(name, to_return, args)
        print s
        exec s
        return eval('%s' % name)


