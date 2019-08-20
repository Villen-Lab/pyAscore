# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from Ascore cimport Ascore 

cdef class PyAscore:
    cdef Ascore c_ascore

    def __cinit__(self, float window_size = 100.):
        self.c_ascore = Ascore(window_size)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def score(self, np.ndarray[double, ndim=1, mode="c"] masses not None,
                    np.ndarray[double, ndim=1, mode="c"] intensities not None):
        cdef int m_len, i_len
        m_len, i_len = masses.shape[0], intensities.shape[0]
