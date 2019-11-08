#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <unordered_map>
#include <numeric>
#include <string>
#include "Ascore.h"
#include "Spectra.h"
#include "ModifiedPeptide.h"
#include "Util.h"

namespace ptmscoring {

    Ascore::Ascore () {
        score_weights = {0.5, 0.75, 1.0, 1.0, 1.0, 1.0, 0.75, 0.5, 0.25, 0.25};
        float weight_sum = std::accumulate(score_weights.begin(), score_weights.end(), 0);
        for (float & val : score_weights){ val /= weight_sum; } 
    }

    Ascore::~Ascore () {}

    void Ascore::resetInternalState () {
        // Clear scoring containers
        peptide_scores_.clear();
        mod_sig_pos_.clear();
        ascores_.clear();

        // Only add scoring distributions if first initialization or depth has changed
        if (scoring_distributions.size() < binned_spectra_ptr->getNTop()){
            for  (size_t depth = scoring_distributions.size() + 1; 
                  depth <= binned_spectra_ptr->getNTop(); 
                  depth++){
                scoring_distributions.push_back( depth / 100. );
            }
        }
    }

    bool Ascore::isUnambiguous () {

        if (modified_peptide_ptr->getNumberOfMods() == modified_peptide_ptr->getNumberModifiable()){
            peptide_scores_.push_back( {} );
            peptide_scores_.front().signature.resize(modified_peptide_ptr->getNumberOfMods(), 1);
            peptide_scores_.front().weighted_score = std::numeric_limits<float>::infinity();
            ascores_.resize(modified_peptide_ptr->getNumberOfMods(), 
                            std::numeric_limits<float>::infinity());
            return true;
        } else {
            return false;
        }

    }

    void Ascore::accumulateCounts () {
        for (char fragment_type : {'b', 'y'}){
            size_t signature_ind = 0;
            // Initialize count stack with a zero count vector
            std::unordered_map<size_t, std::vector<size_t>> count_map; 
            count_map[1] = std::vector<size_t>(binned_spectra_ptr->getNTop());

            // Iterate through signatures and accumulate scores
            std::vector<size_t> cur_count;
            ModifiedPeptide::FragmentGraph graph = modified_peptide_ptr->getFragmentGraph(
                fragment_type, 1
            );
            for (; !graph.isSignatureEnd(); graph.incrSignature()){

                // Copy last signature state
                cur_count = count_map[graph.getFragmentSize()];

                for (; !graph.isFragmentEnd(); graph.incrFragment()){
                    if ( graph.isPositionModified() ) { 
                        count_map[graph.getFragmentSize()] = cur_count; 
                    }
                    float fragment_mz = graph.getFragmentMZ();
                    if ( modified_peptide_ptr->hasMatch(fragment_mz) ){
                        std::tuple<float, size_t> matched_peak = modified_peptide_ptr->getMatch(
                            fragment_mz
                        );
                        cur_count[std::get<1>(matched_peak)]++;
                    }
                }

                // Write signature info
                if ( fragment_type == 'b' ) {
                    peptide_scores_.push_back({ graph.getSignature(), cur_count, {} , -1});
                } else {
                    ScoreContainer & score_ref = *(peptide_scores_.end() - 1 - signature_ind);
                    for (size_t ind = 0; ind < cur_count.size(); ind++) {
                        score_ref.counts[ind] += cur_count[ind];
                    }
                    signature_ind += 1;
                }
            }
        }

        for (ScoreContainer & scores : peptide_scores_) {
            for (std::vector<size_t>::iterator it = scores.counts.begin() + 1;
                 it != scores.counts.end(); it++) {
                 *it += *(it - 1);
            }
        }   
    }

    void Ascore::calculateFullScores () {
        size_t peptide_size = modified_peptide_ptr->getPeptide().size();

        for ( ScoreContainer & score_cont : peptide_scores_ ) {
            for (size_t nions = 1; nions <= score_cont.counts.size(); nions++) {
                score_cont.scores.push_back(
                    std::abs(
                        -10 * scoring_distributions.at(nions - 1).log10_pvalue(
                            score_cont.counts.at(nions - 1), 2 * peptide_size
                        )
                    )
                );
            }
            score_cont.weighted_score = std::inner_product(
                score_weights.begin(), score_weights.end(), score_cont.scores.begin(), 0.
            );
        }
    }

    void Ascore::sortScores() {
        auto score_greater = [](const ScoreContainer &a, const ScoreContainer &b) {
            return a.weighted_score > b.weighted_score;
        };
        std::sort(peptide_scores_.begin(), peptide_scores_.end(), score_greater);
    }

    void Ascore::findModifiedPos () {
        ScoreContainer & best_score = peptide_scores_.front();
        for (size_t ind = 0; ind < best_score.signature.size(); ind++) {
            if (best_score.signature[ind]) {mod_sig_pos_.push_back(ind);}
        }
    }

    std::vector<size_t> Ascore::findDifferences(const ScoreContainer& ref, 
                                                const ScoreContainer& other) {
        std::vector<size_t> differences;
        for (size_t ind = 0; ind < mod_sig_pos_.size(); ind++) {
            if ( ref.signature[mod_sig_pos_[ind]] != other.signature[mod_sig_pos_[ind]] ) {
                differences.push_back(ind);
            }
        }

        return differences;
    } 

    float Ascore::calculateAmbiguity(const ScoreContainer& ref, const ScoreContainer& other){
        // Remove trivial cases
        if (std::abs(ref.weighted_score - other.weighted_score) < 1e-6) {
            return 0.;
        }

        // Find best depth
        float max_score_diff = 0;
        float max_score_depth = 0;
        for (size_t depth = 0; depth < ref.scores.size(); depth++) {
            float diff = ref.scores[depth] - other.scores[depth];
            if (diff > max_score_diff) {
                max_score_diff = diff;
                max_score_depth = depth;
            }
        }
        //std::cout << max_score_diff << " " << max_score_depth << std::endl;

        // Calculate ascore
        std::vector<size_t> ion_counts(2);
        std::vector<size_t> ion_trials(2);
        for (char fragment_type : {'b', 'y'}) {
            std::vector<std::vector<float>> ions = modified_peptide_ptr->getSiteDeterminingIons(
                ref.signature, other.signature, fragment_type, 1
            );

            ion_trials[0] += ions[0].size();
            for (float mz : ions[0]){
                if (modified_peptide_ptr->hasMatch(mz) and
                    std::get<1>(modified_peptide_ptr->getMatch(mz)) <= max_score_depth){
                    ion_counts[0]++;
                }
            }

            ion_trials[1] += ions[1].size();
            for (float mz : ions[1]){
                if (modified_peptide_ptr->hasMatch(mz) and
                    std::get<1>(modified_peptide_ptr->getMatch(mz)) <= max_score_depth){
                    ion_counts[1]++;
                }
            }   
        }
        
        std::vector<float> scores(2);
        for (size_t ind : {0, 1}) {
            scores[ind] = std::abs(
                -10 * scoring_distributions.at(max_score_depth).log10_pvalue(
                          ion_counts[ind], ion_trials[ind]
                )
            );
        }

        return scores[0] - scores[1];
    }

    void Ascore::calculateAscores() {
        findModifiedPos();
        ascores_.resize(mod_sig_pos_.size(), std::numeric_limits<float>::infinity());

        ScoreContainer & best_score = peptide_scores_.front();
        //std::cout << std::endl << "Best:" << std::endl;
        //for ( size_t i : best_score.signature ) {std::cout << i << " ";}
        //std::cout << best_score.weighted_score << std::endl << "Competing:" << std::endl;
        for ( ScoreContainer & competing_score : peptide_scores_ ) {
            const std::vector<size_t> & differences = findDifferences(best_score, competing_score);
            if ( differences.size() == 1 and !std::isfinite(ascores_[differences.front()]) ) {
                ascores_[differences.front()] = calculateAmbiguity(best_score, competing_score);
                //for ( size_t i : competing_score.signature ) {std::cout << i << " ";}
                //std::cout << competing_score.weighted_score << " "
                //<< ascores_[differences.front()] << std::endl;
            }
        }   
    }

    void Ascore::score (const BinnedSpectra & binned_spectra, 
                        const ModifiedPeptide & modified_peptide) {
        // Reinitialize for every peptide that comes in
        modified_peptide_ptr = &modified_peptide;
        binned_spectra_ptr = &binned_spectra;
        resetInternalState();

        // Score pipeline
        if ( not isUnambiguous() ){
            accumulateCounts();
            calculateFullScores();
            sortScores();
            calculateAscores();
        }
    }

    std::string Ascore::getBestSequence () {
        if ( peptide_scores_.size() ) {
            return modified_peptide_ptr->getPeptide(peptide_scores_.front().signature);
        } else {
            return "";
        }
    }

    float Ascore::getBestScore () {
        if ( peptide_scores_.size() ) {
            return peptide_scores_.front().weighted_score;
        } else {
            return -1;
        }
    }

    std::vector<float> Ascore::getAscores () {
        return ascores_;
    }
}

