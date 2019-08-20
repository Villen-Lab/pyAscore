cdef extern from "lib/Spectra.cpp":
    pass

cdef extern from "lib/Spectra.h" namespace "ptmscoring":
    cdef cppclass BinnedSpectra:
        BinnedSpectra(float, float, float, int) except +
        void consumeSpectra(const double *, const double *, int)

        const double & getMZ()
        const double & getIntensity()
        size_t getNPeaks()

        size_t getBin()
        size_t resetBin()
        size_t setBin(size_t)
        size_t nextBin()

        size_t getRank()
        size_t resetRank()
        size_t setRank(size_t)
        size_t nextRank()

        float getMinMZ()
        float getMaxMZ()
        float getBinSize()

        size_t getNBins()
        size_t getNTop()

