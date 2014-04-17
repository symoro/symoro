#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This module performs writing and reading data into PAR file. PAR is a
plain text file used to represent the different parameters of the robot.
"""


import os
import re

from symoroutils import filemgr
from symoroutils import tools
from pysymoro import robot


_keywords = ['ant', 'sigma', 'b', 'd', 'r',
             'gamma', 'alpha', 'mu', 'theta',
             'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ',
             'MX', 'MY', 'MZ', 'M',
             'IA', 'FV', 'FS', 'FX', 'FY', 'FZ',
             'CX', 'CY', 'CZ', 'QP', 'QDP', 'GAM',
             'W0', 'WP0', 'V0', 'VP0',
             'Z', 'G']

_NF = ['ant', 'sigma', 'b', 'd', 'r',
       'gamma', 'alpha', 'mu', 'theta']
_NJ = ['QP', 'QDP', 'GAM']
_NL = ['XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ',
       'MX', 'MY', 'MZ', 'M',
       'IA', 'FV', 'FS', 'FX', 'FY', 'FZ',
       'CX', 'CY', 'CZ']
_VEC = ['W0', 'WP0', 'V0', 'VP0']
_ZERO_BASED = {'W0', 'WP0', 'V0', 'VP0', 'Z', 'G'}
_bool_dict = {'True': True, 'False': False,
              'true': True, 'false': False,
              '1': True, '0': False}

_keyword_repl = {'Ant': 'ant', 'Sigma': 'sigma', 'B': 'b', 'R': 'r',
                 'Alpha': 'alpha', 'Mu': 'mu', 'Theta': 'theta'}


def _extract_vals(robo, key, line):
    line = line.replace('{', '')
    line = line.replace('}', '')
    if key in _ZERO_BASED:
        k = 0
    else:
        k = 1
    items = line.split(',')
    items_proc = []
    prev_item = False
    for i, v in enumerate(items):
        if v.find('atan2') == -1 and not prev_item:
            items_proc.append(v)
        elif prev_item:
            items_proc.append('%s,%s' % (items[i-1], v))
            prev_item = False
        else:
            prev_item = True
    for i, v in enumerate(items_proc):
        if robo.put_val(i+k, key, v.strip()) == tools.FAIL:
            return tools.FAIL


def _write_par_list(robo, f, key, N0, N):
    f.write('%s = {%s' % (key, robo.get_val(N0, key)))
    for i in xrange(N0 + 1, N):
        f.write(',%s' % robo.get_val(i, key))
    f.write('}\n')


def writepar(robo):
    fname = filemgr.make_file_path(robo)
    with open(fname, 'w') as f:
        f.write('(* Robotname = \'%s\' *)\n' % robo.name)
        f.write('NL = %s\n' % robo.nl)
        f.write('NJ = %s\n' % robo.nj)
        f.write('NF = %s\n' % robo.nf)
        f.write('Type = %s\n' % tools.TYPES.index(robo.structure))
        f.write('is_mobile = %s\n' % int(robo.is_mobile))
        f.write('\n(* Geometric parameters *)\n')
        if robo.is_mobile:
            N0 = 0
        else:
            N0 = 1
        for key in _NF:
            _write_par_list(robo, f, key, 1, robo.NF)
        f.write('\n(* Dynamic parameters and external forces *)\n')
        for key in _NL:
            _write_par_list(robo, f, key, N0, robo.NL)
        f.write('\n(* Joint parameters *)\n')
        for key in _NJ:
            _write_par_list(robo, f, key, 1, robo.NJ)
        f.write('\n(* Speed and acceleration of the base *)\n')
        for key in _VEC:
            _write_par_list(robo, f, key, 0, 3)
        f.write('\n(* Acceleration of gravity *)\n')
        _write_par_list(robo, f, 'G', 0, 3)
        f.write('\n(* Transformation of 0 frame position fT0 *)\n')
        _write_par_list(robo, f, 'Z', 0, 16)
        f.write('\n(* End of definition *)\n')


def readpar(robo_name, file_path):
    """Return:
        robo: an instance of Robot, read from file
        flag: indicates if any errors accured. (tools.FAIL)
    """
    with open(file_path, 'r') as f:
        #initialize the Robot instance
        f.seek(0)
        d = {}
        is_mobile = False
        for line in f.readlines():
            for s in ('NJ', 'NL', 'Type'):
                match = re.match(r'^%s.*=([\d\s]*)(\(\*.*)?' % s, line)
                if match:
                    d[s] = int(match.group(1))
                    continue
            match = re.match(r'^is_mobile.*=([\d\s]*)(\(\*.*)?', line)
            if match:
                is_mobile = _bool_dict[(match.group(1).strip())]
        if len(d) < 2:
            return None, tools.FAIL
        NF = d['NJ']*2 - d['NL']
        robo = robot.Robot(robo_name, d['NL'], d['NJ'], NF,
                            is_mobile, tools.TYPES[d['Type']])
        robo.directory = os.path.dirname(file_path)
        #fitting the data
        acc_line = ''
        key = ''
        f.seek(0)
        flag = tools.OK
        for line in f.readlines():
            if line.find('(*') != -1:
                continue
            line = line.replace('Pi', 'pi')
            match = re.match(r'^(.*)=.*\{(.*)', line)
            if match:
                acc_line == ''
                key = match.group(1).strip()
                acc_line = match.group(2).strip()
            else:
                acc_line += line
            if acc_line.find('}') != -1:
                if key in _keyword_repl:
                    key = _keyword_repl[key]
                if key in _keywords:
                    if _extract_vals(robo, key, acc_line) == tools.FAIL:
                        flag = tools.FAIL
                acc_line = ''
                key = ''
    return robo, flag


