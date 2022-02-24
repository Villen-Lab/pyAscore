from setuptools import setup, Extension
from Cython.Build import cythonize
from os import path

import numpy as np

NAME = "pyascore"
VERSION = 0.2
DESCR = "A python module for fast post translational modification localization."
REQUIRES = ['cython']

AUTHOR = "Anthony Valente"
EMAIL = "valenta4@uw.edu"

LICENSE = "MIT"

SRC_DIR = "pyascore"
PACKAGES = [SRC_DIR, SRC_DIR + ".parsing"]

HERE = path.abspath(path.dirname(__file__))
# Get the long description from the README file
with open(path.join(HERE, 'README.md'), encoding='utf-8') as f:
    LONG_DESCR = f.read()

EXT = [Extension(
          SRC_DIR + ".ptm_scoring", 
          [SRC_DIR + "/ptm_scoring/PTMScoring.pyx"],
          include_dirs=[np.get_include(), "/ptm_scoring/lib"],
          extra_compile_args=["/std:c++latest"]
      )]

setup(
    name=NAME,
    version=VERSION,
    description=DESCR,
    long_description=LONG_DESCR,
    long_description_content_type='text/markdown',
    packages=PACKAGES,
    install_requires=REQUIRES,
    python_requires='>=3.5',
    ext_modules=cythonize(EXT, compiler_directives={'language_level' : "3"}),
    entry_points={
          'console_scripts': [
              'pyascore = pyascore.__main__:main'
          ]
    }
) 
