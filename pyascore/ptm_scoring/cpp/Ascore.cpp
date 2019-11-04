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
        best_weighted_score = -1;

        // Clear scoring containers
        signatures.clear();
        accumulated_counts.clear();
        peptide_scores.clear();

        // Clear distributions and reinitialize
        scoring_distributions.clear();
        for (size_t depth = 1; depth <= binned_spectra_ptr->getNTop(); depth++){
            scoring_distributions.push_back( depth / 100. );
        }
    }

    bool Ascore::isUnambiguous () {

        if (modified_peptide_ptr->getNumberOfMods() == modified_peptide_ptr->getNumberModifiable()){
            best_signature = std::vector<size_t>(modified_peptide_ptr->getNumberOfMods());
            std::fill(best_signature.begin(), best_signature.end(), 1);

            best_weighted_score = std::numeric_limits<float>::infinity();
            return true;
        } else {
            return false;
        }

    }

    void Ascore::accumulateCounts () {
        for (char fragment_type : {'b'}){
            size_t signature_ind = 0;
            // Initialize count stack with a zero count vector
            std::unordered_map<size_t, std::vector<size_t>> count_map; 
            count_map[1] = std::vector<size_t>(binned_spectra_ptr->getNTop());

            // Iterate through signatures and accumulate scores
            std::vector<size_t> cur_count;
            ModifiedPeptide::FragmentGraph graph = modified_peptide_ptr->getFragmentGraph(fragment_type, 1);
            for (;
                 !graph.isSignatureEnd();
                 graph.incrSignature()){
                //for (size_t i : graph.getSignature() ){
                //    std::cout << i << " ";
                //}
                //std::cout << std::endl << graph.getFragmentSeq() << std::endl << std::endl;

                // Copy last signature state
                cur_count = count_map[graph.getFragmentSize()];

                for (; !graph.isFragmentEnd(); graph.incrFragment()){
                    if ( graph.isPositionModified() ) { 
                        count_map[graph.getFragmentSize()] = cur_count; 
                    }
                    float fragment_mz = graph.getFragmentMZ();
                    if ( modified_peptide_ptr->hasMatch(fragment_mz) ){
                        std::tuple<float, size_t> matched_peak = modified_peptide_ptr->getMatch(fragment_mz);
                        cur_count[std::get<1>(matched_peak)]++;
                    }
                }

                // Write signature info
                if ( fragment_type == 'b' ) {
                    signatures.push_back(graph.getSignature());
                    accumulated_counts.push_back(cur_count);
                } else {
                    std::vector<size_t> & count_ref = *(accumulated_counts.end() - 1 - signature_ind);
                    for (size_t ind = 0; ind < cur_count.size(); ind++) {
                        count_ref[ind] += cur_count[ind];
                    }
                    signature_ind += 1;
                }
            }
        }

        for (std::vector<size_t> & counts : accumulated_counts) {
            for (std::vector<size_t>::iterator it = counts.begin() + 1;
                 it != counts.end(); it++) {
                 *it += *(it - 1);
            }
        }   
    }

    void Ascore::calculateFullScores () {
        size_t peptide_size = modified_peptide_ptr->getPeptide().size();
        //std::cout << modified_peptide_ptr->getPeptide() << std::endl;

        std::vector<std::vector<size_t>>::iterator sig_it = signatures.begin();
        std::vector<std::vector<size_t>>::iterator count_it = accumulated_counts.begin();
        for (; sig_it != signatures.end(); sig_it++, count_it++) {
            peptide_scores.push_back({});
            for (size_t nions = 1; nions <= count_it->size(); nions++) {
                //std::cout << count_it->at(nions - 1) << " ";
                peptide_scores.back().push_back(
                    std::abs(
                        -10 * scoring_distributions.at(nions - 1).log10_pvalue(
                            count_it->at(nions - 1), 2 * peptide_size
                        )
                    )
                );       
            }
            peptide_scores.back().push_back(
                std::inner_product(score_weights.begin(), score_weights.end(), peptide_scores.back().begin(), 0.)
            );
            //std::cout << peptide_scores.back().back() << " ";
            //std::cout << std::endl;
        }
    }

    void Ascore::findBestPeptide () {
        size_t best_score_ind = 0;
        best_weighted_score = peptide_scores.front().back();
        for (auto score_it = peptide_scores.begin(); score_it != peptide_scores.end(); score_it++) {
            if (score_it->back() > best_weighted_score){
                best_weighted_score = score_it->back();
                best_score_ind = score_it - peptide_scores.begin();
            }
        }
        best_signature = signatures[best_score_ind];
    }

    void Ascore::score (const BinnedSpectra & binned_spectra, const ModifiedPeptide & modified_peptide) {
        // Reinitialize for every peptide that comes in
        modified_peptide_ptr = &modified_peptide;
        binned_spectra_ptr = &binned_spectra;
        resetInternalState();

        // Score pipeline
        if ( not isUnambiguous() ){
            accumulateCounts();
            calculateFullScores();
            findBestPeptide();
        }
    }

    std::string Ascore::getBestSequence () {
        if ( best_weighted_score >= 0 ) {
            return modified_peptide_ptr->getPeptide(best_signature);
        } else {
            return "";
        }
    }

    float Ascore::getBestScore () {
        return best_weighted_score;
    }

}

