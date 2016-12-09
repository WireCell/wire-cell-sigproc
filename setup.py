#!/usr/bin/env python
'''
Build wirecell.sigproc package.

FIXME: this currently is broken if installed in the wider context of a
top-level wirecell package.
'''

from setuptools import setup, find_packages
setup(
    name = 'wirecell-sigproc',
    provides = ["wirecell.sigproc"],
    version = '0.0',
    package_dir = {'':'python'},
    packages = ['wirecell', 'wirecell.sigproc', 'wirecell.sigproc.response', 'wirecell.sigproc.paper'],
    install_requires = [
        'Click',
    ],
    entry_points = dict(
        console_scripts = [
            'wirecell-sigproc = wirecell.sigproc.main:main',
        ]
    )
)

