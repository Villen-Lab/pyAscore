#ifndef ASCORE_H
#define ASCORE_H

#include <vector>
#include <stack>
#include <string>
#include "Spectra.h"
#include "ModifiedPeptide.h"
#include "Util.h"

namespace ptmscoring {

    struct ScoreContainer {
        std::vector<size_t> signature;
        std::vector<size_t> counts;
        std::vector<float> scores;
        float weighted_score = -1;
        size_t total_fragments;
    };

    struct AscoreContainer {
        size_t sig_pos; // Position within signature of reference
        std::vector<size_t> competing_index; // Position within peptide of competing modifiable amino acids
        std::vector<float> pep_scores; // Pepscore for competing modifiable amino acids
        std::vector<float> ascores; // Ascores of reference vs other modifiable amino acids
    };

    class Ascore {
        BinomialDist bin_dist;
        const BinnedSpectra * binned_spectra_ptr;
        const ModifiedPeptide * modified_peptide_ptr;

        std::vector<float> score_weights;
        std::vector<BinomialDist> scoring_distributions;

        std::vector<ScoreContainer> peptide_scores_;
        std::vector<AscoreContainer> ascore_containers_;

        void resetInternalState();
        void accumulateCounts();
        void calculateFullScores();
        void sortScores();
        bool isUnambiguous();
        void findModifiedPos();
        void calculateAscores();
        public:
            Ascore();
            ~Ascore();

            void score(const BinnedSpectra &, const ModifiedPeptide &);
            float calculateAmbiguity(const ScoreContainer&, const ScoreContainer&);
            std::string getBestSequence();
            std::vector<std::string> getAllSequences();
            float getBestScore();
            std::vector<ScoreContainer> getAllPepScores();
            std::vector<float> getAscores();
            std::vector<size_t> getAlternativeSites(size_t);
    };

}

#endif
