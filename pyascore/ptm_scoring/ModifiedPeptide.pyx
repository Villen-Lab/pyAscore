# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
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

    def reset_iterator(self, char fragment_type):
        self.modified_peptide_ptr[0].resetIterator(fragment_type)

    def incr_signature(self):
        return self.modified_peptide_ptr[0].incrSignature()

    def get_signature(self):
        cdef vector[size_t] signature_vector = self.modified_peptide_ptr[0].getSignature()
        signature_array = np.zeros(signature_vector.size(), dtype=np.uint64)

        cdef size_t i = 0
        for i in range(signature_vector.size()):
            signature_array[i] = signature_vector[i]
        return signature_array

    def incr_fragment(self):
        return self.modified_peptide_ptr[0].incrFragment()

    def get_fragment_type(self):
        return self.modified_peptide_ptr[0].getFragmentType()

    def get_fragment_mz(self, size_t charge):
        return self.modified_peptide_ptr[0].getFragmentMZ(charge)

    def get_fragment_size(self):
        return self.modified_peptide_ptr[0].getFragmentSize()

    def get_fragment_seq(self):
        return self.modified_peptide_ptr[0].getFragmentSeq().decode('utf8')
