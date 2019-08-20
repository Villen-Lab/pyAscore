# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from libcpp.string cimport string
from ModifiedPeptide cimport ModifiedPeptide

cdef class PyModifiedPeptide:
    cdef ModifiedPeptide * modified_peptide_ptr

    def __cinit__(self, str mod_group, float mod_mass):
        self.modified_peptide_ptr = new ModifiedPeptide(mod_group.encode("utf8"), mod_mass)

    def __dealloc__(self):
        del self.modified_peptide_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def consume_peptide(self, str peptide, size_t n_of_mod, 
                        np.ndarray[size_t, ndim=1, mode="c"] aux_mod_pos = None, 
                        np.ndarray[float, ndim=1, mode="c"] aux_mod_mass = None):

        self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod)
