#!/usr/bin/env python
# -*- coding: utf-8 -*-


from setuptools import setup, find_packages


def readme():
    with open('README.md') as f:
        return f.read()


setup(
    name='symoro',
    version='0.1alpha',
    description='SYmoblic MOdelling of RObots software',
    url='http://github.com/vijaravind/symoro',
    scripts=['bin/symoro-bin'],
    packages=find_packages(),
    zip_safe=False
)


