# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from Ascore cimport Ascore
from ModifiedPeptide cimport ModifiedPeptide
from Spectra cimport BinnedSpectra

cdef class PyAscore:
    cdef Ascore * ascore_ptr
    cdef ModifiedPeptide * modified_peptide_ptr
    cdef BinnedSpectra * binned_spectra_ptr

    def __cinit__(self, float min_mz, float max_mz, 
                        float bin_size, size_t n_top,
                        str mod_group, float mod_mass):
        self.binned_spectra_ptr = new BinnedSpectra(min_mz, max_mz, bin_size, n_top)
        self.modified_peptide_ptr = new ModifiedPeptide(mod_group.encode("utf8"), mod_mass)
        self.ascore_ptr = new Ascore()

    def __dealloc__(self):
        del self.binned_spectra_ptr
        del self.modified_peptide_ptr
        del self.ascore_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def score(self, np.ndarray[double, ndim=1, mode="c"] mz_arr not None, 
                    np.ndarray[double, ndim=1, mode="c"] int_arr not None,
                    str peptide, size_t n_of_mod):
        cdef size_t m_len, i_len
        m_len, i_len = mz_arr.size, int_arr.size
        self.binned_spectra_ptr[0].consumeSpectra(&mz_arr[0], &int_arr[0], m_len)

        self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod)
        while (self.binned_spectra_ptr[0].getBin() < self.binned_spectra_ptr[0].getNBins()):

            self.binned_spectra_ptr[0].resetRank()
            while (self.binned_spectra_ptr[0].getRank() < self.binned_spectra_ptr[0].getNPeaks()):
                self.modified_peptide_ptr[0].consumePeak(self.binned_spectra_ptr[0].getMZ(),
                                                         self.binned_spectra_ptr[0].getRank())
                self.binned_spectra_ptr[0].nextRank()

            self.binned_spectra_ptr[0].nextBin()

        self.ascore_ptr[0].score(self.binned_spectra_ptr[0], self.modified_peptide_ptr[0])
