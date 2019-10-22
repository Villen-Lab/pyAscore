#ifndef ASCORE_H
#define ASCORE_H

#include <vector>
#include <stack>
#include "Spectra.h"
#include "ModifiedPeptide.h"
#include "Util.h"

namespace ptmscoring {

    class Ascore {
        BinomialDist bin_dist;
        const BinnedSpectra * binned_spectra_ptr;
        const ModifiedPeptide * modified_peptide_ptr;

        std::vector<std::vector<size_t>> signatures;
        std::vector<std::vector<size_t>> accumulated_counts;
        std::vector<std::vector<float>> peptide_scores;
        std::vector<float> score_weights;
        std::vector<BinomialDist> scoring_distributions;

        std::stack<std::vector<size_t>> count_stack;
        void accumulate_counts();
        void calculate_full_scores();
        public:
            Ascore();
            ~Ascore();

            void score(const BinnedSpectra &, const ModifiedPeptide &);
    };

}

#endif
