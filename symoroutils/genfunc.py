# -*- coding: utf-8 -*-
"""
This module provides MATLAB function generation
"""

from sympy import Matrix, Symbol


def gen_fheader_matlab(symo, name, args,
                       multival=False):
    func_head = []
    func_head.append('function RESULT=%s (' % name)
    func_head.append(convert_syms_matlab(args))
    func_head.append(')\n')
    return func_head


def convert_mat_matlab(to_return):
    sym_list = convert_to_list(to_return, keep_const=True)
    res = []
    sym_iter = iter(sym_list)
    for i in range(to_return.shape[0]):
        for j in range(to_return.shape[1]):
            res.append(str(sym_iter.next()))
            res.append(',')
        res[-1] = ';'
    res.pop()
    return "".join(res)


def convert_syms_matlab(syms):
    """Converts 'syms' structure to a string

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
            res.extend(convert_to_list(item, keep_const))
        return res
    elif isinstance(syms, Symbol) or keep_const:
        return [syms]
    else:
        return []


def gen_fbody_matlab(symo, name, to_return, args, ret_name=''):
    """Generates list of string statements of the function that
    computes symbolf from to_return.  arg_syms are considered to
    be known
    """
     # set of defined symbols
    arg_syms = symo.extract_syms(args)
    # final symbols to be compute
    multiline_res = False   # may be useful for C/C++ code
    if len(to_return) > 16 and isinstance(to_return, Matrix):
        multiline_res = True
        to_return_list = [symo.sydi[s] for s in to_return if s in symo.sydi]
        res_syms = symo.extract_syms(to_return_list)
    else:
        res_syms = symo.extract_syms(to_return)
        if isinstance(to_return, Matrix):
            to_ret_str = convert_mat_matlab(to_return)
        else:
            to_ret_str = convert_syms_matlab(to_return)
    # defines order of computation
    order_list = symo.sift_syms(res_syms, arg_syms)
    # list of instructions in final function
    func_body = []
    # will be switched to true when branching detected
    space = '    '
    folded = 1    # indentation = 1 + number of 'for' statements
    multival = False
    glob = 0
    glob_item = ''
    for s in order_list:
        if s not in symo.sydi:
            if glob == 0:
                glob += 1
                glob_item += '%sglobal %s' % (space * folded, s)
            elif glob < 12:
                glob_item += ' %s' % s
                glob += 1
            else:
                glob = 0
                glob_item += ' %s\n' % s
        else:
            if isinstance(symo.sydi[s], tuple):
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
            item = item.replace('**', '^')
            func_body.append(item)
    if multiline_res:
        rows, cols = to_return.shape
        func_body.insert(0, '%sRESULT=zeros(%s,%s);\n' % (space, rows, cols))
        form_str = space + 'RESULT(%s,%s)=%s;\n'
        for i in range(rows):
            for j in range(cols):
                s = to_return[i, j]
                if s in symo.sydi:
                    item = form_str % (i + 1, j + 1, symo.sydi[s])
                    item = item.replace('**', '^')
                    func_body.append(item)
    elif multival:
        func_body.insert(0, '%sRESULT=[];\n' % (space))
        item = '%sRESULT=[RESULT;%s];\n' % (space * folded, to_ret_str)
        func_body.append(item)
    else:
        item = '%sRESULT=[%s];\n' % (space * folded, to_ret_str)
        func_body.append(item)
    for f in range(folded-1, 0, -1):
        func_body.append('%send\n' % (space * f))
    func_body.append('end\n')
    func_body.insert(0, glob_item + '\n')
    return func_body
