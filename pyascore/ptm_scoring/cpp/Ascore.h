#ifndef ASCORE_H
#define ASCORE_H

#include <vector>
#include <stack>
#include <string>
#include "Spectra.h"
#include "ModifiedPeptide.h"
#include "Util.h"

namespace ptmscoring {

    class Ascore {
        BinomialDist bin_dist;
        const BinnedSpectra * binned_spectra_ptr;
        const ModifiedPeptide * modified_peptide_ptr;

        std::vector<float> score_weights;
        std::vector<BinomialDist> scoring_distributions;

        struct ScoreContainer {
            std::vector<size_t> signature;
            std::vector<size_t> counts;
            std::vector<float> scores;
            float weighted_score = -1;
        };
        std::vector<ScoreContainer> peptide_scores_;

        std::vector<size_t> mod_sig_pos_;
        std::vector<float> ascores_;

        void resetInternalState();
        void accumulateCounts();
        void calculateFullScores();
        void sortScores();
        bool isUnambiguous();
        void findModifiedPos();
        std::vector<size_t> findDifferences(const ScoreContainer&, const ScoreContainer&);
        float calculateAmbiguity(const ScoreContainer&, const ScoreContainer&);
        void calculateAscores();
        public:
            Ascore();
            ~Ascore();

            void score(const BinnedSpectra &, const ModifiedPeptide &);
            std::string getBestSequence();
            float getBestScore();
            std::vector<float> getAscores();
    };

}

#endif
