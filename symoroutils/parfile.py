# -*- coding: utf-8 -*-


# This file is part of the OpenSYMORO project. Please see
# https://github.com/symoro/symoro/blob/master/LICENCE for the licence.


"""
This module performs writing and reading data into PAR file. PAR is a
plain text file used to represent the different parameters of the robot.
"""


import os
import re

from symoroutils import filemgr
from symoroutils import tools
from pysymoro import robot


_keywords = [
    'ant', 'sigma', 'b', 'd', 'r', 'gamma', 'alpha', 'mu', 'theta',
    'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ', 'MX', 'MY', 'MZ', 'M',
    'IA', 'FV', 'FS', 'FX', 'FY', 'FZ', 'CX', 'CY', 'CZ',
    'eta', 'k', 'QP', 'QDP', 'GAM', 'W0', 'WP0', 'V0', 'VP0', 'Z', 'G'
]
_NF = ['ant', 'sigma', 'b', 'd', 'r', 'gamma', 'alpha', 'mu', 'theta']
_NJ = ['eta', 'k', 'QP', 'QDP', 'GAM']
_NL = [
    'XX', 'XY', 'XZ', 'YY', 'YZ', 'ZZ', 'MX', 'MY', 'MZ', 'M',
    'IA', 'FV', 'FS', 'FX', 'FY', 'FZ', 'CX', 'CY', 'CZ'
]
_VEC = ['W0', 'WP0', 'V0', 'VP0']
_ZERO_BASED = {'W0', 'WP0', 'V0', 'VP0', 'Z', 'G'}
_bool_dict = {
    'True': True,
    'False': False,
    'true': True,
    'false': False,
    '1': True,
    '0': False
}
_keyword_repl = {
    'Ant': 'ant',
    'Mu': 'mu',
    'Sigma': 'sigma',
    'B': 'b',
    'Alpha': 'alpha',
    'Theta': 'theta',
    'R': 'r'
}


def _extract_vals(robo, key, line):
    line = line.replace('{', '')
    line = line.replace('}', '')
    if key in _ZERO_BASED:
        k = 0
    elif (robo.is_floating or robo.is_mobile) and key in _NL:
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
    f.write('{0} = {{{1}'.format(key, robo.get_val(N0, key)))
    for i in range(N0 + 1, N):
        f.write(',{0}'.format(robo.get_val(i, key)))
    f.write('}\n')


def writepar(robo):
    fname = robo.par_file_path
    with open(fname, 'w') as f:
        # robot description
        f.write('(* Robotname = \'{0}\' *)\n'.format(robo.name))
        f.write('NL = {0}\n'.format(robo.nl))
        f.write('NJ = {0}\n'.format(robo.nj))
        f.write('NF = {0}\n'.format(robo.nf))
        f.write('Type = {0}\n'.format(tools.TYPES.index(robo.structure)))
        f.write('is_floating = {0}\n'.format(robo.is_floating))
        f.write('is_mobile = {0}\n'.format(robo.is_mobile))
        # geometric parameters
        f.write('\n(* Geometric parameters *)\n')
        for key in _NF:
            _write_par_list(robo, f, key, 1, robo.NF)
        # dynamic parameters
        f.write('\n(* Dynamic parameters and external forces *)\n')
        N0 = 0 if robo.is_floating or robo.is_mobile else 1
        for key in _NL:
            _write_par_list(robo, f, key, N0, robo.NL)
        # joint parameters
        f.write('\n(* Joint parameters *)\n')
        for key in _NJ:
            _write_par_list(robo, f, key, 1, robo.NJ)
        # base parameters - velocity and acceleration
        f.write('\n(* Velocity and acceleration of the base *)\n')
        for key in _VEC:
            _write_par_list(robo, f, key, 0, 3)
        # gravity vector
        f.write('\n(* Acceleration of gravity *)\n')
        _write_par_list(robo, f, 'G', 0, 3)
        # base parameters - Z matrix
        f.write('\n(* Transformation of 0 frame position fT0 *)\n')
        _write_par_list(robo, f, 'Z', 0, 16)
        f.write('\n(* End of definition *)\n')


def readpar(robo_name, file_path):
    """Return:
        robo: an instance of Robot, read from file
        flag: indicates if any errors accured. (tools.FAIL)
    """
    with open(file_path, 'r') as f:
        f.seek(0)
        d = {}
        is_floating = False
        is_mobile = False
        for line in f.readlines():
            # check for robot name
            name_pattern = r"\(\*.*Robotname.*=.*\'([\s\w-]*)\'.*\*\)"
            match = re.match(name_pattern, line)
            if match:
                robo_name = match.group(1).strip()
            # check for joint numbers, link numbers, type
            for s in ('NJ', 'NL', 'Type'):
                match = re.match(r'^%s.*=([\d\s]*)(\(\*.*)?' % s, line)
                if match:
                    d[s] = int(match.group(1))
                    continue
            # check for is_floating
            match = re.match(r'^is_floating.*=([\w\s]*)(\(\*.*)?', line)
            if match:
                is_floating = _bool_dict[(match.group(1).strip())]
            # check for is_mobile
            match = re.match(r'^is_mobile.*=([\w\s]*)(\(\*.*)?', line)
            if match:
                is_mobile = _bool_dict[(match.group(1).strip())]
        if len(d) < 2:
            return None, tools.FAIL
        NF = d['NJ']*2 - d['NL']
        #initialize the Robot instance
        robo = robot.Robot(
            name=robo_name,
            NL=d['NL'], NJ=d['NJ'], NF=NF,
            structure=tools.TYPES[d['Type']],
            is_floating=is_floating,
            is_mobile=is_mobile,
            directory=os.path.dirname(file_path),
            par_file_path=file_path
        )
        # fitting the data
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


