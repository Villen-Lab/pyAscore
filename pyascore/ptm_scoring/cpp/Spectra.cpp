#include <iostream>
#include <cmath>
#include <algorithm>
#include <vector>
#include "Spectra.h"
#include "Types.h"

namespace ptmscoring {

    BinnedSpectra::BinnedSpectra (float bin_size, size_t n_top) {
        this->bin_size = bin_size;
        this->n_top = n_top;

        this->min_mz = 0.;
        this->max_mz = std::numeric_limits<float>::infinity();
        this->n_bins = -1;

        resetBin();
        resetRank();
    }

    BinnedSpectra::~BinnedSpectra () {}

    void BinnedSpectra::sortTopSpectra() {
        auto peak_greater = [](const Peak &a, const Peak &b) { return a.intensity > b.intensity; };

        for (size_t bin_ind = 0; bin_ind < spec_bins.size(); bin_ind++) {
            if (n_top < spec_bins[bin_ind].size()){
                // Gather the n_top in unsorted order to the front of the vector
                std::nth_element(spec_bins[bin_ind].begin(), 
                                 spec_bins[bin_ind].begin() + n_top - 1, 
                                 spec_bins[bin_ind].end(), peak_greater);

                // Resize the final bin vector
                spec_bins[bin_ind].resize(n_top);
            } 

            // Sort the top peaks
            std::sort(spec_bins[bin_ind].begin(), spec_bins[bin_ind].end(), peak_greater);
        }
    }

    void BinnedSpectra::consumeSpectra (const double * mz_arr, const double * int_arr, size_t n_peaks) {

        // Decide on bounds and reset internal state
        min_mz = std::floor( *std::min_element(mz_arr, mz_arr + n_peaks) / 100.) * 100.;
        max_mz = std::ceil( *std::max_element(mz_arr, mz_arr + n_peaks) / 100.) * 100.;
        n_bins = std::ceil( (max_mz - min_mz) / bin_size );
        spec_bins.clear();
        spec_bins.resize(n_bins);

        // Copy peaks to bins
        int target_bin;
        for (size_t ind = 0; ind < n_peaks; ind++){
            target_bin = std::min( 
                (size_t) std::floor( ( mz_arr[ind] - min_mz ) / bin_size ), 
                n_bins-1
            );
            spec_bins[target_bin].push_back({mz_arr[ind], int_arr[ind]});
        }

        // Find and maintain only top peaks
        sortTopSpectra();

        // Reset iterator
        resetBin();
        resetRank();
    }

    // This class is setup so that you if you accidently
    // look at a peak beyond the end of a given bin
    // then you get an out_of_range exception.
    const double & BinnedSpectra::getMZ() const {
        return spec_bins.at(bin).at(rank).mz;
    }
    
    const double & BinnedSpectra::getIntensity() const {
        return spec_bins.at(bin).at(rank).intensity;
    }

    size_t BinnedSpectra::getNPeaks() const {
        return spec_bins.at(bin).size();
    }

    // Internal bin state getters and setters
    size_t BinnedSpectra::getBin () const { return bin; }
    size_t BinnedSpectra::resetBin () { return setBin(0); }
    size_t BinnedSpectra::setBin (size_t new_bin) {
        bin = std::min(new_bin, n_bins); 
        return bin;
    }
    size_t BinnedSpectra::nextBin () {
        bin = std::min(bin + 1, n_bins);
        return bin;
    }

    // Internal rank state getters and setters
    size_t BinnedSpectra::getRank () const { return rank; }
    size_t BinnedSpectra::resetRank () { return setRank(0); }
    size_t BinnedSpectra::setRank (size_t new_rank) {
        rank = std::min(new_rank, n_top); 
        return rank;
    }
    size_t BinnedSpectra::nextRank () {
        rank = std::min(rank + 1, n_top);
        return rank;
    }

    // Getters of internal properties
    float BinnedSpectra::getMinMZ () const { return min_mz; }
    float BinnedSpectra::getMaxMZ () const { return max_mz; }
    float BinnedSpectra::getBinSize () const { return bin_size; }

    size_t BinnedSpectra::getNBins () const { return n_bins; }
    size_t BinnedSpectra::getNTop () const { return n_top; }

}
