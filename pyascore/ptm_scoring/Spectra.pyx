# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from Spectra cimport BinnedSpectra

cdef class PyBinnedSpectra:
    """
    The PyBinnedSpectra object shuttles MS/MS peaks to equal sized mass bins

    Bin size is determined at object creation, but lower and upper limits of the bin range
    are set as the highest and lowest MZ of the consumed spectra. After allocating each peak to
    its respective bin, the peaks are sorted and at most n_top of the most intense peakse are retained.
    The object then offers iteration over bins and top peaks.

    Parameters
    ----------
    bin_size : float
        Size in MZ of each bin
    n_top : int
        Number of top peaks to retain in each bin (must be >= 0)

    Attributes
    ----------
    bin_size : int
        Bin size in MZ
    min_mz : float
        Minimum MZ of consumed spectra
    max_mz : float
        Maximum MZ of consumed spectra
    n_bins : int
        Total number of bins able to be fit between min_mz and max_mz 
    mz : float
        Current peak's MZ
    intensity : float
        Current peaks's intensity
    n_peaks : int
        Current bin's number of peaks retained (may be <= ntop)
    bin : int
        Current bin index
    rank : int
        Current peak's rank in bin
    """
    cdef BinnedSpectra * binned_spectra_ptr

    def __cinit__(self, float bin_size, size_t n_top):
        self.binned_spectra_ptr = new BinnedSpectra(bin_size, n_top)

    def __dealloc__(self):
        del self.binned_spectra_ptr

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def consume_spectra(self, np.ndarray[double, ndim=1, mode="c"] mz_arr not None,
                              np.ndarray[double, ndim=1, mode="c"] int_arr not None):
        """
        Consumes the MZ and intensities of a single spectra

        After consumption, bin and rank will be initialized to 0.

        Parameters
        ----------
        mz_arr : ndarray of float64
            Array of MZ values for each peak in a spectra
        int_arr : ndarray of float64
            Array of intensity values for each peak in a spectra
        """
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
        """Sets bin to 0"""
        self.binned_spectra_ptr[0].resetBin()
    def next_bin(self):
        """Increments bin by 1"""
        self.binned_spectra_ptr[0].nextBin()

    @property
    def rank(self):
        return self.binned_spectra_ptr[0].getRank()
    @rank.setter
    def rank(self, new_rank):
        self.binned_spectra_ptr[0].setRank(new_rank)
    def reset_rank(self):
        """Sets rank to 0"""
        self.binned_spectra_ptr[0].resetRank()
    def next_rank(self):
        """Increments rank by 1"""
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
