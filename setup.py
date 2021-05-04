#!/usr/bin/env python
# -*- coding: utf-8 -*-
from setuptools import setup, Extension
from setuptools.command.build_ext import build_ext
from Cython.Build import cythonize
import subprocess
import numpy

class build_ext_make(build_ext):
    def run(self):
        subprocess.check_call('make')
        build_ext.run(self)


optical_depth = Extension(
    'optical_depth',
    sources=['optical_depth.pyx'],
    libraries=['optical_depth', 'gfortran', 'm'],
    library_dirs=['build'],
    include_dirs=[numpy.get_include()],
)

setup(
    name='optical_depth',
    version='1.0',
    author='Alexander Sandrock',
    author_email='alexander.sandrock@tu-dortmund.de',
    ext_modules=cythonize([optical_depth]),
    cmdclass={
        'build_ext': build_ext_make,
    }
)

