#ifndef SPECTRA_H
#define SPECTRA_H

#include <vector>
#include "Types.h"

namespace ptmscoring {
    class BinnedSpectra { 
        private:
            std::vector<std::vector<Peak>> spec_bins;
            float min_mz, max_mz, bin_size;
            size_t n_bins, n_peaks, n_top;
            size_t bin, rank;
            void sortTopSpectra();
        public:
            BinnedSpectra(float, size_t);
            ~BinnedSpectra();
            void consumeSpectra(const double *, const double *, size_t);

            const double & getMZ() const;
            const double & getIntensity() const;
            size_t getNPeaks() const;

            size_t getBin() const;
            size_t resetBin();
            size_t setBin(size_t);
            size_t nextBin();

            size_t getRank() const;
            size_t resetRank();
            size_t setRank(size_t);
            size_t nextRank();

            float getMinMZ() const;
            float getMaxMZ() const;
            float getBinSize() const;

            size_t getNBins() const;
            size_t getNTop() const;

    };
}

#endif
