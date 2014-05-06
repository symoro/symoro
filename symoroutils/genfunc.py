# -*- coding: utf-8 -*-
"""
This module provides MATLAB function generation
"""

from sympy import Matrix, Symbol
from sympy import sympify


def gen_fheader_matlab(symo, name, args, to_return, multival=False):
    func_head = ['function [']
    if multival:
        func_head.append(convert_syms_matlab(to_return, 'r'))
    else:
        func_head.append(convert_syms_matlab(to_return))
    func_head.append('] = %s_func (' % name)
    func_head.append(convert_syms_matlab(args))
    func_head.append(')\n');
    return func_head

def convert_syms_matlab(syms, extension=''):
    """Converts 'syms' structure to sintactically correct string

    Parameters
    ==========
    syms: list, Matrix or tuple of them
    rpl_liter: bool
        if true, all literals will be replaced with _
        It is done to evoid expression like [x, 0] = args[1]
        Because it will cause exception of assigning to literal
    """
#    if isinstance(syms, tuple) or isinstance(syms, list):
#        syms = [convert_syms_matlab(item, rpl_liter) for item in syms]
#        res = []
#        for i, s in enumerate(syms):
#            res.append(s)
#            if i < len(syms) - 1:
#                res += ','
#        return res
#    elif isinstance(syms, Matrix):
#        res = ''
#        for i in xrange(syms.shape[0]):
#            res += convert_syms_matlab(list(syms[i, :]), rpl_liter)
#            if i < syms.shape[0] - 1:
#                res += ','
#        return res
#    elif rpl_liter and sympify(syms).is_number:
#        return ''
#    else:
#        return '%s%s' % (syms, extension)
    cond1 = isinstance(syms, tuple)
    cond2 = isinstance(syms, list)
    cond3 = isinstance(syms, Matrix)
    if cond1 or cond2 or cond3:
        res = []
        for item in iter(syms):
            s = convert_syms_matlab(item, extension)
            if s is not None:
                res.append(s)
                res.append(',')
        res.pop()
        return "".join(res)
    elif isinstance(syms, Symbol):
        return '%s%s' % (syms, extension)
    else:
        return None

def convert_to_list(syms):
    cond1 = isinstance(syms, tuple)
    cond2 = isinstance(syms, list)
    cond3 = isinstance(syms, Matrix)
    if cond1 or cond2 or cond3:
        res = []
        for item in syms:
            res.extend(convert_to_list(item))
        return res
    elif isinstance(syms, Symbol):
            return syms
    else:
        return []


def gen_fbody_matlab(symo, name, to_return, args):
    """Generates list of string statements of the function that
    computes symbolf from to_return.  wr_syms are considered to
    be known
    """
     # set of defined symbols
    wr_syms = symo.extract_syms(args)
    # final symbols to be compute
    syms = symo.extract_syms(to_return)
    syms_list = list(syms)
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
        item.replace('**','^')
        func_body.append(item)
    if multival:
        for sym in syms_list:
            func_body.insert(0, '%s%sr=[];\n' % (space, sym))
            item = '%s%sr=[%sr;%s];\n' % (space*folded, sym, sym, sym)
            func_body.append(item)

    for f in xrange(folded-1, 0, -1):
        func_body.append('%send\n' % (space * f))
    func_body.append('end\n')

    return func_body, multival
