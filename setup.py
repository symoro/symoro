#!/usr/bin/env python
# -*- coding: utf-8 -*-


import os
from setuptools import setup, find_packages


BIN_FOLDER = 'bin'


def readme():
    with open('README.md') as f:
        return f.read()


def apply_folder_join(item):
    return os.path.join(BIN_FOLDER, item)


if os.name is 'nt':
    bin_scripts = ['symoro-bin.py']
else:
    bin_scripts = ['symoro-bin']
bin_scripts = map(apply_folder_join, bin_scripts)


setup(
    name='symoro',
    version='0.2',
    description='SYmoblic MOdelling of RObots software package',
    url='http://github.com/symoro/symoro',
    license='MIT',
    scripts=bin_scripts,
    packages=find_packages(exclude=['*.tests', '*.tests.*', 'tests.*', 'tests']),
    install_requires=[
        'sympy>=0.7.3',
        'numpy>=1.6.1',
        'wxPython>=2.8.11',
        'PyOpenGL>=3.0.1b2'
    ],
    zip_safe=False
)


