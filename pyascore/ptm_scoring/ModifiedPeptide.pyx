# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.string cimport string
from ModifiedPeptide cimport ModifiedPeptide

cdef class PyModifiedPeptide:
    cdef ModifiedPeptide * modified_peptide_ptr

    def __cinit__(self, str mod_group, float mod_mass, float mz_error = .5):
        self.modified_peptide_ptr = new ModifiedPeptide(mod_group.encode("utf8"), mod_mass, mz_error)

    def __dealloc__(self):
        del self.modified_peptide_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def consume_peptide(self, str peptide, size_t n_of_mod, 
                        np.ndarray[unsigned int, ndim=1, mode="c"] aux_mod_pos = None, 
                        np.ndarray[float, ndim=1, mode="c"] aux_mod_mass = None):
        if aux_mod_pos is not None and aux_mod_mass is not None:
            self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod,
                                                        &aux_mod_pos[0], &aux_mod_mass[0], 
                                                        aux_mod_pos.size)
        else:
            self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod)

    def get_fragment_graph(self, str fragment_type, size_t charge_state):
        return PyFragmentGraph(self, fragment_type.encode("utf8")[0], charge_state)

    def get_site_determining_ions(self, np.ndarray[unsigned int, ndim=1, mode="c"] sig_1,
                                        np.ndarray[unsigned int, ndim=1, mode="c"] sig_2,
                                        str fragment_type, size_t charge_state):

        cdef vector[size_t] sig_vec_1
        cdef vector[size_t] sig_vec_2
        cdef size_t ind
        for ind in range(<size_t> min(sig_1.size, sig_2.size)):
            sig_vec_1.push_back(<size_t> sig_1[ind])
            sig_vec_2.push_back(<size_t> sig_2[ind])

        cdef vector[vector[float]] ions = self.modified_peptide_ptr[0].getSiteDeterminingIons(
            sig_vec_1, sig_vec_2, fragment_type.encode("utf8")[0], charge_state
        )

        ion_arrays = (np.zeros(ions[0].size(), dtype=np.float32),
                      np.zeros(ions[1].size(), dtype=np.float32))
        for ind in range(ions[0].size()):
            ion_arrays[0][ind] = ions[0][ind]
        for ind in range(ions[1].size()):
            ion_arrays[1][ind] = ions[1][ind]

        return ion_arrays

cdef class PyFragmentGraph:
    cdef ModifiedPeptide.FragmentGraph * fragment_graph_ptr

    def __cinit__(self, PyModifiedPeptide peptide, char fragment_type, size_t charge_state):
        self.fragment_graph_ptr = new ModifiedPeptide.FragmentGraph(
            peptide.modified_peptide_ptr, fragment_type, charge_state
        )

    def __dealloc__(self):
        del self.fragment_graph_ptr

    @property
    def fragment_type(self):
        return self.fragment_graph_ptr[0].getFragmentType()

    @property
    def charge_state(self):
        return self.fragment_graph_ptr[0].getChargeState()

    def reset_iterator(self):
        return self.fragment_graph_ptr[0].resetIterator()

    def incr_signature(self):
        return self.fragment_graph_ptr[0].incrSignature()

    def is_signature_end(self):
        return self.fragment_graph_ptr[0].isSignatureEnd()

    def incr_fragment(self):
        return self.fragment_graph_ptr[0].incrFragment()

    def is_fragment_end(self):
        return self.fragment_graph_ptr[0].isFragmentEnd()

    def set_signature(self, np.ndarray[unsigned int, ndim=1, mode="c"] new_signature):
        cdef vector[size_t] signature_vector
        cdef size_t ind
        for ind in range(<size_t> new_signature.size):
            signature_vector.push_back(<size_t> new_signature[ind])

        self.fragment_graph_ptr[0].setSignature(signature_vector)
        

    def get_signature(self):
        cdef vector[size_t] signature_vector = self.fragment_graph_ptr[0].getSignature()
        signature_array = np.zeros(signature_vector.size(), dtype=np.uint64)

        cdef size_t i = 0
        for i in range(signature_vector.size()):
            signature_array[i] = signature_vector[i]
        return signature_array

    def get_fragment_mz(self):
        return self.fragment_graph_ptr[0].getFragmentMZ()

    def get_fragment_size(self):
        return self.fragment_graph_ptr[0].getFragmentSize()

    def get_fragment_seq(self):
        return self.fragment_graph_ptr[0].getFragmentSeq().decode('utf8')

