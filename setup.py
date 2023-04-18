from setuptools import setup, Extension
from Cython.Build import cythonize
from os import path, environ
from sys import platform

import numpy as np

extra_compile_args=[]
if "PYASCORE_COMPILE_ARGS" in environ:
    extra_compile_args=environ["PYASCORE_COMPILE_ARGS"].split()

if platform == "darwin":
    extra_compile_args.append('-std=c++11')

SRC_DIR = "pyascore"
EXT = [
    Extension(SRC_DIR + ".ptm_scoring",
              [SRC_DIR + "/ptm_scoring/PTMScoring.pyx"],
              include_dirs=[np.get_include()],
              extra_compile_args=extra_compile_args)
]

setup(ext_modules=cythonize(EXT, compiler_directives={'language_level' : "3"}))
