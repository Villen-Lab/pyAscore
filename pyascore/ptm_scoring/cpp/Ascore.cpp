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
        float weight_sum = std::accumulate(score_weights.begin(), score_weights.end(), 0.);
        for (float & val : score_weights){ val /= weight_sum; } 
    }

    Ascore::~Ascore () {}

    void Ascore::resetInternalState () {
        // Clear scoring containers
        peptide_scores_.clear();
        ascore_containers_.clear();

        // Only add scoring distributions if first initialization or depth has changed
        if (scoring_distributions.size() < binned_spectra_ptr->getNTop()){
            for  (size_t depth = scoring_distributions.size() + 1; 
                  depth <= binned_spectra_ptr->getNTop(); 
                  depth++){
                scoring_distributions.push_back( 2 * modified_peptide_ptr->getMZError() * depth / 100. );
            }
        }
    }

    bool Ascore::isUnambiguous () {

        if (modified_peptide_ptr->getNumberOfMods() >= modified_peptide_ptr->getNumberModifiable()){
            for (size_t ind = 0; ind < modified_peptide_ptr->getNumberOfMods(); ind++){
                ascore_containers_.push_back({
                    ind, {}, {}, {std::numeric_limits<float>::infinity()}
                });
            }
            return true;
        } else {
            return false;
        }

    }

    void Ascore::accumulateCounts () {
        std::unordered_map<long, ScoreContainer> score_cache;
        for (char fragment_type : modified_peptide_ptr->getFragmentTypes()){
	    for (size_t charge = 1; charge <= modified_peptide_ptr->getMaxFragmentCharge(); charge++) {
                // Initialize count stack with a zero count vector
                std::unordered_map<size_t, std::vector<size_t>> count_cache;
                std::unordered_map<size_t, size_t> nfrag_cache;
                count_cache[1] = std::vector<size_t>(binned_spectra_ptr->getNTop());
                nfrag_cache[1] = 0;

                // Iterate through signatures and accumulate scores
                std::vector<size_t> cur_count;
                size_t cur_nfrag;
                ModifiedPeptide::FragmentGraph graph = modified_peptide_ptr->getFragmentGraph(
                    fragment_type, charge
                );
                for (; !graph.isSignatureEnd(); graph.incrSignature()){

                    // Copy last signature state
                    cur_count = count_cache[graph.getFragmentSize()];
                    cur_nfrag = nfrag_cache[graph.getFragmentSize()];

                    for (; !graph.isFragmentEnd(); graph.incrFragment()){
                        if ( graph.isPositionModified() && !graph.isLoss() ) { 
                            count_cache[graph.getFragmentSize()] = cur_count;
                            nfrag_cache[graph.getFragmentSize()] = cur_nfrag;
                        }
                        float fragment_mz = graph.getFragmentMZ();
                        if ( modified_peptide_ptr->hasMatch(fragment_mz) ){
                            std::tuple<float, size_t> matched_peak = modified_peptide_ptr->getMatch(
                                fragment_mz
                            );
                            cur_count[std::get<1>(matched_peak)]++;
                        }
                        cur_nfrag++;
                    }

                    // Write signature info
                    long index = 0;
                    for (size_t sig_pos : graph.getSignature()) {
                        index = (index<<1) | sig_pos;
                    }

                    if (!score_cache.count(index)) {
                        score_cache[index] = {};
                        score_cache[index].signature = graph.getSignature();
                        score_cache[index].counts = cur_count;
                        score_cache[index].total_fragments = cur_nfrag;
                        //score_cache[index] = { graph.getSignature(), cur_count, {} , -1, cur_nfrag};
                    } else {
                        ScoreContainer & score_ref = score_cache[index];
                        for (size_t ind = 0; ind < cur_count.size(); ind++) {
                            score_ref.counts[ind] += cur_count[ind];
                        }
                        score_ref.total_fragments += cur_nfrag;
                    }
                }
            }
	}
        
        peptide_scores_.reserve(score_cache.size());
        for (auto & cached_pair : score_cache) {
            for (std::vector<size_t>::iterator it = cached_pair.second.counts.begin() + 1;
                 it != cached_pair.second.counts.end(); it++) {
                 *it += *(it - 1);
            }
            peptide_scores_.push_back(std::move(cached_pair.second));
        }
    }

    void Ascore::calculateFullScores () {
        //size_t peptide_size = modified_peptide_ptr->getBasePeptide().size();
        for ( ScoreContainer & score_cont : peptide_scores_ ) {
            for (size_t nions = 1; nions <= score_cont.counts.size(); nions++) {
                score_cont.scores.push_back(
                    std::abs(
                        -10 * scoring_distributions.at(nions - 1).log10_pvalue(
                            score_cont.counts.at(nions - 1), score_cont.total_fragments
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
            if (best_score.signature[ind]) {
                ascore_containers_.push_back({ind});
            }
        }
    }

    float Ascore::calculateAmbiguity(const ScoreContainer& ref, const ScoreContainer& other){
        // Remove trivial cases
        if (std::abs(ref.weighted_score - other.weighted_score) < 1e-6) {
            return 0.;
        }

        // Find best depth
        float max_score_diff = 0.;
        size_t max_score_depth = 0;
        for (size_t depth = 0; depth < ref.scores.size(); depth++) {
            float diff = ref.scores[depth] - other.scores[depth];
            if (diff > max_score_diff) {
                max_score_diff = diff;
                max_score_depth = depth;
            }
        }

        // Calculate ascore
        std::vector<size_t> ion_counts(2);
        std::vector<size_t> ion_trials(2);
        for (char fragment_type : modified_peptide_ptr->getFragmentTypes()) {
            std::vector<std::vector<float>> ions = modified_peptide_ptr->getSiteDeterminingIons(
                ref.signature, other.signature, fragment_type,
		modified_peptide_ptr->getMaxFragmentCharge()
            );

            ion_trials[0] += ions[0].size();
            for (float mz : ions[0]){
                if (modified_peptide_ptr->hasMatch(mz) &&
                    std::get<1>(modified_peptide_ptr->getMatch(mz)) <= max_score_depth){
                    ion_counts[0]++;
                }
            }

            ion_trials[1] += ions[1].size();
            for (float mz : ions[1]){
                if (modified_peptide_ptr->hasMatch(mz) &&
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

        ScoreContainer & best_score = peptide_scores_.front();
        for ( ScoreContainer & competing_score : peptide_scores_ ) {
            size_t ndifferences = modified_peptide_ptr->getNumberOfMods() - std::inner_product(
                best_score.signature.begin(), best_score.signature.end(), 
                competing_score.signature.begin(), 0
            );

            if ( ndifferences != 1 ) continue;

            size_t same_count = 0;
            size_t ascore_ind = 0;
            size_t comp_sig_ind = 0;
            std::vector<size_t>::iterator best_it = best_score.signature.begin();
            std::vector<size_t>::iterator comp_it = competing_score.signature.begin();
            for (; best_it < best_score.signature.end(); best_it++, comp_it++) {
                if (*best_it && *comp_it) {
                    same_count++;
                } else if (*best_it && (!*comp_it)) {
                    ascore_ind = same_count;
                } else if ((! *best_it) && *comp_it) {
                    comp_sig_ind = comp_it - competing_score.signature.begin();
                }
            }

            AscoreContainer & cont = ascore_containers_.at(ascore_ind);
            if (cont.ascores.empty() || 
                competing_score.weighted_score == cont.pep_scores.back() ) {
                cont.competing_index.push_back(
                    modified_peptide_ptr->getPosOfNthModifiable(comp_sig_ind) + 1 
                );
                cont.pep_scores.push_back( 
                    competing_score.weighted_score 
                );
                cont.ascores.push_back( 
                    calculateAmbiguity(best_score, competing_score) 
                );
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
	accumulateCounts();
        calculateFullScores();
	    
        if ( !isUnambiguous() ){
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

    std::vector<std::string> Ascore::getAllSequences () {
        std::vector<std::string> sequences;
        for (ScoreContainer cont : peptide_scores_) {
            sequences.push_back(modified_peptide_ptr->getPeptide(cont.signature));
        }
        return sequences;
    }

    float Ascore::getBestScore () {
        if ( peptide_scores_.size() ) {
            return peptide_scores_.front().weighted_score;
        } else {
            return -1;
        }
    }

    std::vector<ScoreContainer> Ascore::getAllPepScores () {
        if ( peptide_scores_.size() ) {
            return peptide_scores_;
        } else {
            return {};
        }
    }

    std::vector<float> Ascore::getAscores () {
        std::vector<float> ascores; 
        for (AscoreContainer& cont : ascore_containers_) {
            ascores.push_back(
                *std::min_element(cont.ascores.begin(), cont.ascores.end())
            );
        }
        return ascores;
    }

    std::vector<size_t> Ascore::getAlternativeSites (size_t site) {
        std::vector<size_t> ambiguity_range = ascore_containers_.at(site).competing_index;
        std::sort(ambiguity_range.begin(), ambiguity_range.end());
        return ambiguity_range;
    }

}

