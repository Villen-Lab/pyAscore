# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from Spectra cimport BinnedSpectra

cdef class PyBinnedSpectra:
    cdef BinnedSpectra * binned_spectra_ptr

    def __cinit__(self, float bin_size, size_t n_top):
        self.binned_spectra_ptr = new BinnedSpectra(bin_size, n_top)

    def __dealloc__(self):
        del self.binned_spectra_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def consume_spectra(self, np.ndarray[double, ndim=1, mode="c"] mz_arr not None,
                              np.ndarray[double, ndim=1, mode="c"] int_arr not None):
        cdef size_t m_len, i_len
        m_len, i_len = mz_arr.size, int_arr.size
        self.binned_spectra_ptr[0].consumeSpectra(&mz_arr[0], &int_arr[0], m_len)

    # Peak access
    @property
    def mz(self):
        return self.binned_spectra_ptr[0].getMZ()

    @property
    def intensity(self):
        return self.binned_spectra_ptr[0].getIntensity()

    @property
    def n_peaks(self):
        return self.binned_spectra_ptr[0].getNPeaks()

    # Mutable state access
    @property
    def bin(self):
        return self.binned_spectra_ptr[0].getBin()
    @bin.setter
    def bin(self, new_bin):
        self.binned_spectra_ptr[0].setBin(new_bin)
    def reset_bin(self):
        self.binned_spectra_ptr[0].resetBin()
    def next_bin(self):
        self.binned_spectra_ptr[0].nextBin()

    @property
    def rank(self):
        return self.binned_spectra_ptr[0].getRank()
    @rank.setter
    def rank(self, new_rank):
        self.binned_spectra_ptr[0].setRank(new_rank)
    def reset_rank(self):
        self.binned_spectra_ptr[0].resetRank()
    def next_rank(self):
        self.binned_spectra_ptr[0].nextRank()

    # Immutable State access
    @property
    def min_mz(self):
        return self.binned_spectra_ptr[0].getMinMZ()
    @property
    def max_mz(self):
        return self.binned_spectra_ptr[0].getMaxMZ()
    @property
    def bin_size(self):
        return self.binned_spectra_ptr[0].getBinSize()
    @property
    def n_bins(self):
        return self.binned_spectra_ptr[0].getNBins()
