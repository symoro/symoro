# -*- coding: utf-8 -*-
"""
This module provides MATLAB function generation
"""

from sympy import Matrix, Symbol


def gen_fheader_matlab(symo, name, args,
                       multival=False):
    func_head = []
    func_head.append('function RESULT=%s_func (' % name)
    func_head.append(convert_syms_matlab(args))
    func_head.append(')\n')
    return func_head


def convert_mat_matlab(to_return):
    sym_list = convert_to_list(to_return, keep_const=True)
    res = []
    sym_iter = iter(sym_list)
    for i in xrange(to_return.shape[0]):
        for j in xrange(to_return.shape[1]):
            res.append(str(sym_iter.next()))
            res.append(',')
        res[-1] = ';'
    res.pop()
    return "".join(res)


def convert_syms_matlab(syms):
    """Converts 'syms' structure to sintactically correct string

    Parameters
    ==========
    syms: list, Matrix or tuple of them
    rpl_liter: bool
        if true, all literals will be replaced with _
        It is done to evoid expression like [x, 0] = args[1]
        Because it will cause exception of assigning to literal
    """
    sym_list = convert_to_list(syms, keep_const=False)
    res = []
    for item in iter(sym_list):
        res.append(str(item))
        res.append(',')
    res.pop()
    return "".join(res)


def convert_to_list(syms, keep_const=True):
    cond1 = isinstance(syms, tuple)
    cond2 = isinstance(syms, list)
    cond3 = isinstance(syms, Matrix)
    if cond1 or cond2 or cond3:
        res = []
        for item in syms:
            res.extend(convert_to_list(item))
        return res
    elif isinstance(syms, Symbol) or keep_const:
        return [syms]
    else:
        return []


def gen_fbody_matlab(symo, name, to_return, args, ret_name=''):
    """Generates list of string statements of the function that
    computes symbolf from to_return.  wr_syms are considered to
    be known
    """
     # set of defined symbols
    wr_syms = symo.extract_syms(args)
    # final symbols to be compute
    syms = symo.extract_syms(to_return)
    if isinstance(to_return, Matrix):
        to_ret_str = convert_mat_matlab(to_return)
    else:
        to_ret_str = convert_syms_matlab(to_return)
    # defines order of computation
    order_list = symo.sift_syms(syms, wr_syms)
    # list of instructions in final function
    func_body = []
    # will be switched to true when branching detected
    space = '    '
    folded = 1    # indentation = 1 + number of 'for' statements
    multival = False
    for s in order_list:
        if s not in symo.sydi:
            item = '%s%s=1.;\n' % (space * folded, s)
        elif isinstance(symo.sydi[s], tuple):
            multival = True
            items = ['%sfor %s=' % (space * folded, s)]
            for x in symo.sydi[s]:
                items.append('%s' % x)
                items.append(',')
            items.append('\n')
            item = "".join(items)
            folded += 1
        else:
            item = '%s%s=%s;\n' % (space * folded, s, symo.sydi[s])
        item.replace('**', '^')
        func_body.append(item)
    if multival:
        func_body.insert(0, '%sRESULT=[];\n' % (space))
        item = '%sRESULT=[RESULT;%s];\n' % (space * folded, to_ret_str)
        func_body.append(item)
    else:
        item = '%sRESULT=[%s];\n' % (space * folded, to_ret_str)
        func_body.append(item)
    for f in xrange(folded-1, 0, -1):
        func_body.append('%send\n' % (space * f))
    func_body.append('end\n')

    return func_body
