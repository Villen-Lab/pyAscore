#include <iostream>
#include <cmath>
#include <vector>
#include <tuple>
#include <stack>
#include <numeric>
#include "Ascore.h"
#include "Spectra.h"
#include "ModifiedPeptide.h"
#include "Util.h"

namespace ptmscoring {

    Ascore::Ascore () {
        score_weights = {0.5, 0.75, 1.0, 1.0, 1.0, 
                         1.0, 0.75, 0.5, 0.25, 0.25}; 
    }

    Ascore::~Ascore () {}

    void Ascore::accumulate_counts () {
        for (char fragment_type : {'b', 'y'}){
            size_t signature_ind = 0;
            // Initialize count stack with a zero count vector
            count_stack = std::stack<std::vector<size_t>>(); 
            count_stack.emplace( std::vector<size_t>(binned_spectra_ptr->getNTop()) );

            // Iterate through signatures and accumulate scores
            std::vector<size_t> cur_count;
            for (ModifiedPeptide::FragmentGraph graph = modified_peptide_ptr->getFragmentGraph(fragment_type, 1);
                 !graph.isSignatureEnd(); 
                 graph.incrSignature()){

                // Copy last signature state
                cur_count = count_stack.top();
                count_stack.pop();

                for (; !graph.isFragmentEnd(); graph.incrFragment()){
                    if ( graph.isPositionModified() ) { count_stack.push(cur_count); }
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


    void Ascore::calculate_full_scores () {
        size_t peptide_size = modified_peptide_ptr->getPeptide().size();
        std::vector<std::vector<size_t>>::iterator sig_it = signatures.begin();
        std::vector<std::vector<size_t>>::iterator count_it = accumulated_counts.begin();
        for (; sig_it != signatures.end(); sig_it++, count_it++) {
            peptide_scores.push_back({});
            for (size_t nions = 1; nions <= count_it->size(); nions++) {
                peptide_scores.back().push_back(
                    -10 * scoring_distributions.at(nions - 1).log_pvalue(
                        count_it->at(nions - 1), 2 * peptide_size
                    )
                );
            }
            peptide_scores.back().push_back(
                std::inner_product(score_weights.begin(), score_weights.end(), peptide_scores.back().begin(), 0.)
            );
        }
    }

    void Ascore::score (const BinnedSpectra & binned_spectra, const ModifiedPeptide & modified_peptide) {
        std::cout << std::endl << "Scoring peptide" << std::endl;

        modified_peptide_ptr = &modified_peptide;
        binned_spectra_ptr = &binned_spectra;

        accumulate_counts();

        for (size_t depth = 1; depth <= binned_spectra_ptr->getNTop(); depth++){
            scoring_distributions.push_back( depth / 100. );
        }
        calculate_full_scores();


        std::cout << "Finished peptide" << std::endl;
    }
}

