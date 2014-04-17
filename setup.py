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
    version='0.1alpha',
    description='SYmoblic MOdelling of RObots software',
    url='http://github.com/symoro/symoro',
    scripts=bin_scripts,
    packages=find_packages(),
    zip_safe=False
)


