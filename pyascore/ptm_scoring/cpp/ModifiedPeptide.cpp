#include <iostream>
#include <cstdio>
#include <vector>
#include <tuple>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <ctype.h>
#include "Types.h"
#include "ModifiedPeptide.h"

namespace ptmscoring {

    ModifiedPeptide::ModifiedPeptide (std::string mod_group, float mod_mass, float mz_error, std::string fragment_types){
        this->mod_group = mod_group;
        this->mod_mass = mod_mass;
        this->mz_error = mz_error;
        this->fragment_types = fragment_types;
    }

    ModifiedPeptide::~ModifiedPeptide () {};

    void ModifiedPeptide::initializeResidues () {

        // Dealloc if already full
        residues.clear();
        neutral_losses.clear();

        // Iterate through characters of peptide
        for (std::string::iterator res_it = peptide.begin();
	     res_it != peptide.end();
	     res_it ++) {
            residues.push_back({ StandardResidues.at(*res_it) });

            // If character is in neutral loss group, add neutral loss
            neutral_losses.push_back({0.});
            if (nl_groups.count(*res_it)) {
                neutral_losses.back().front() = nl_groups.at(*res_it);
            }

            // If character is in mod group, add a modified residue
	    // Also, if terminal mods are allowed and is terminal aa, add a modified residue
	    bool is_mod_aa = mod_group.find(*res_it) != std::string::npos;
	    bool allow_nterm_mod = mod_group.find('n') != std::string::npos && res_it == peptide.begin();
	    bool allow_cterm_mod = mod_group.find('c') != std::string::npos && res_it + 1 == peptide.end();
            if (is_mod_aa || allow_nterm_mod || allow_cterm_mod) {
                residues.back().push_back( residues.back().back() + mod_mass );
                neutral_losses.back().push_back(0.);
                if (nl_groups.count(std::tolower(*res_it))) {
                    neutral_losses.back().back() = nl_groups.at(std::tolower(*res_it));
                } 
            }

        }

    }

    void ModifiedPeptide::applyAuxMods () {
        for (size_t ind = 0; ind < aux_mod_pos.size(); ind++) {
	    // Correct indices
	    size_t pep_ind = aux_mod_pos[ind];
	    if (pep_ind > 0) {
                pep_ind -= 1;
	    }

            // Apply aux mod to both modified and unmodified
	    for (float& m : residues[pep_ind]) {
                m += aux_mod_mass[ind];
            }

	    // Add NL groups
	    char res = peptide.at(pep_ind);
            if (nl_groups.count(std::tolower(res))) {
                neutral_losses[pep_ind].front() = nl_groups.at(std::tolower(res));
                neutral_losses[pep_ind].resize(1);
            }
	}
    }

    void ModifiedPeptide::initializeFragments () {
        fragment_scores.clear();
        fragments.clear();
        for (char t : fragment_types){
            for (size_t charge=1; charge <= max_fragment_charge; charge++) {
                for (FragmentGraph graph = getFragmentGraph(t, charge);
                     !graph.isSignatureEnd();
                     graph.incrSignature()) {
                    for(; !graph.isFragmentEnd();
                          graph.incrFragment()) {
                        fragments.push_back(graph.getFragmentMZ());
                    }
                }
	    }
        }
        std::sort(fragments.begin(), fragments.end());
    }

    void ModifiedPeptide::addNeutralLoss(std::string group, float mass) {\
        for (char const& aa : group) {
           nl_groups[aa] = mass;
        }
    }

    void ModifiedPeptide::consumePeptide (std::string peptide,
		                          size_t n_of_mod,
					  size_t max_fragment_charge,
                                          const unsigned int * aux_mod_pos,
                                          const float * aux_mod_mass,
                                          size_t n_aux_mods) {

        this->peptide = peptide;
        this->n_of_mod = n_of_mod;
	this->max_fragment_charge = max_fragment_charge;

        this->aux_mod_pos.resize(n_aux_mods);
        std::copy(aux_mod_pos, aux_mod_pos + n_aux_mods, this->aux_mod_pos.begin());
        this->aux_mod_mass.resize(n_aux_mods);
        std::copy(aux_mod_mass, aux_mod_mass + n_aux_mods, this->aux_mod_mass.begin());

        initializeResidues();
        applyAuxMods();
        initializeFragments();
    }

    void ModifiedPeptide::consumePeak (float mz, size_t rank) {
        std::vector<float>::iterator it;
        it = lower_bound(fragments.begin(), fragments.end(), mz - .5);

        float ppm_low;
        float ppm_high;
        for(; it < fragments.end(); it++) {
            ppm_low = *it - mz_error;
            ppm_high = *it + mz_error;
            if ((mz > ppm_low) && (mz < ppm_high)) {
		    if ((!fragment_scores.count(*it)) || (std::get<1>(fragment_scores.at(*it)) > rank)) {
			    fragment_scores[*it] = {mz, rank}; // The theoretical mz needs to be the key
		    }
            } else if (mz < ppm_low) {break;}
        }

    }

    bool ModifiedPeptide::hasMatch (float mz) const {
        return fragment_scores.count(mz);
    }

    std::tuple<float, size_t> ModifiedPeptide::getMatch (float mz) const {
        return fragment_scores.at(mz);
    }
    
    std::string ModifiedPeptide::getModGroup () const {
        return mod_group;
    }

    float ModifiedPeptide::getModMass () const { 
        return mod_mass;
    }

    float ModifiedPeptide::getMZError () const {
        return mz_error;
    }

    std::string ModifiedPeptide::getFragmentTypes () const {
        return fragment_types;
    }

    size_t ModifiedPeptide::getNumberOfMods () const {
        return n_of_mod;
    }

    size_t ModifiedPeptide::getMaxFragmentCharge () const {
        return max_fragment_charge;
    }

    size_t ModifiedPeptide::getNumberModifiable () const {
        size_t n_modifiable = 0;
        for (const std::vector<float> & res_it : residues) {
            if (res_it.size() > 1) {n_modifiable += 1;} 
        }
        return n_modifiable;
    }

    size_t ModifiedPeptide::getPosOfNthModifiable (size_t n) const {
        for (size_t ind = 0; ind < residues.size(); ind++) {
            if ((residues.at(ind).size() > 1) && (n == 0)) {
                return ind;
            } else if (residues.at(ind).size() > 1) {
                n--;
            }
        }
        return residues.size();
    }

    std::string ModifiedPeptide::getBasePeptide () const {
        return peptide;
    }

    std::string ModifiedPeptide::getPeptide(std::vector<size_t> signature) const {
        // If no signature is given, just use first signature
        if ( signature.size() == 0 ) {
	    size_t fake_sig_size = std::min(n_of_mod, getNumberModifiable());
            signature = std::vector<size_t>(fake_sig_size, 1);
        }

        std::vector<float> all_mod_masses = std::vector<float>(peptide.size() + 2, 0.);

	// Add overflow mod
	if (n_of_mod > getNumberModifiable()) {
	    if (mod_group.find("n") != std::string::npos) {
                all_mod_masses.front() += mod_mass;
	    } else {
                all_mod_masses.back() +=  mod_mass;
	    }
	}

        // Add internal variable mods
	for (size_t sig_ind = 0; sig_ind < signature.size(); sig_ind++) {
            if (signature[sig_ind] == 1) {
	        size_t mod_pos = getPosOfNthModifiable(sig_ind);
		char aa = peptide[mod_pos];
		if (mod_group.find(aa) != std::string::npos) {
                    all_mod_masses[mod_pos + 1] += mod_mass;
		} else if (mod_pos == 0) {
		    all_mod_masses.front() += mod_mass;
		} else if (mod_pos + 1 == peptide.size()) {
                    all_mod_masses.back() += mod_mass;
		}
	    }
	}

	// Add internal fixed mods
	for (size_t aux_ind = 0; aux_ind < aux_mod_pos.size(); aux_ind++) {
            all_mod_masses[aux_mod_pos[aux_ind]] += aux_mod_mass[aux_ind];
	}

	// Write peptide
	size_t start_ind = 0;
	size_t end_ind = all_mod_masses.size();
	if (all_mod_masses.front() == 0.) start_ind++;
	if (all_mod_masses.back() == 0.) end_ind--;
	std::string full_peptide = "n" + peptide + "c";     
        std::string mod_peptide = "";
	for (size_t ind = start_ind; ind < end_ind; ind++) {
	    mod_peptide += full_peptide[ind];
	    if (all_mod_masses[ind] > 0.) {
                char mod_buffer[10];
                std::sprintf(mod_buffer, "[%d]", (int) std::round(all_mod_masses[ind]));
                mod_peptide += mod_buffer;
	    }
	}
        return mod_peptide;
    }

    ModifiedPeptide::FragmentGraph ModifiedPeptide::getFragmentGraph (char fragment_type, size_t charge_state) const {
        return ModifiedPeptide::FragmentGraph(this, fragment_type, charge_state);
    }

    std::vector<std::vector<float>> ModifiedPeptide::getSiteDeterminingIons (const std::vector<size_t> & signature_1,
                                                                             const std::vector<size_t> & signature_2,
                                                                             char fragment_type, size_t max_charge) const {
	// Generate candidate fragments for each charge state up to max
        std::vector<float> candidate_fragments_1;
	std::vector<float> candidate_fragments_2;
	for (size_t charge = 1; charge <= max_charge; charge++) {
	    // Generate fragments for signature 1
            ModifiedPeptide::FragmentGraph graph_1 = getFragmentGraph(fragment_type, charge);
	    graph_1.setSignature(signature_1);
            for(; !graph_1.isFragmentEnd();
                  graph_1.incrFragment()) {
              candidate_fragments_1.push_back( graph_1.getFragmentMZ() );
	    }

	    // Generate fragments for signature 2
	    ModifiedPeptide::FragmentGraph graph_2 = getFragmentGraph(fragment_type, charge);
            graph_2.setSignature(signature_2);
            for(; !graph_2.isFragmentEnd();
                  graph_2.incrFragment()) {
              candidate_fragments_2.push_back( graph_2.getFragmentMZ() );
            }
        }

        // Sort fragments
	std::sort(candidate_fragments_1.begin(), candidate_fragments_1.end());
        std::sort(candidate_fragments_2.begin(), candidate_fragments_2.end());

        // Determine non-overlapping fragments
	std::vector<std::vector<float>> ion_lists(2);
	std::vector<float>::iterator fragment_iter_1 = candidate_fragments_1.begin();
	std::vector<float>::iterator fragment_iter_2 = candidate_fragments_2.begin();
	while (fragment_iter_1 < candidate_fragments_1.end() ||
	       fragment_iter_2 < candidate_fragments_2.end()) {
	    // 1st and 2nd: 
	    //   Check if one list is just finished. Just grab the rest of the other list.
	    // 3rd:
	    //   Check if fragments are equal given mz_error. These we will discard.
	    // 4th and 5th:
	    //   Check which list to take a fragment from. These are site determining peaks. 
            if (fragment_iter_2 == candidate_fragments_2.end()) {
                ion_lists[0].push_back(*fragment_iter_1);
                fragment_iter_1++;
            } else if (fragment_iter_1 == candidate_fragments_1.end()) {
                ion_lists[1].push_back(*fragment_iter_2);
                fragment_iter_2++;
            } else if (std::abs(*fragment_iter_1 - *fragment_iter_2) < mz_error) {
                fragment_iter_1++;
                fragment_iter_2++;
	    } else if (*fragment_iter_1 < *fragment_iter_2) { 
	        ion_lists[0].push_back(*fragment_iter_1);
		fragment_iter_1++;
	    } else {
	        ion_lists[1].push_back(*fragment_iter_2);
		fragment_iter_2++;
	    }

	}

        return ion_lists;

    }

    ////////////////////////////////////////
    // Start FragmentGraph Implementation //
    ////////////////////////////////////////

    ModifiedPeptide::FragmentGraph::FragmentGraph (const ModifiedPeptide * modified_peptide, char fragment_type, size_t charge_state) {
        this->modified_peptide = modified_peptide;
        this->fragment_type = fragment_type;
        this->charge_state = charge_state;
        resetIterator();
    }

    ModifiedPeptide::FragmentGraph::~FragmentGraph () {}

    char ModifiedPeptide::FragmentGraph::getFragmentType () { return fragment_type; }

    size_t ModifiedPeptide::FragmentGraph::getChargeState () { return charge_state; }

    void ModifiedPeptide::FragmentGraph::resetResidueInd () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            residue_ind = 0;
        } else if ( fragment_type == 'y'  || fragment_type == 'z' || fragment_type == 'Z' ) {
            residue_ind = modified_peptide->residues.size() - 1;
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::setResidueInd (size_t new_ind) {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            residue_ind = std::min(residue_ind, new_ind);
        } else if ( fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z' ) {
            residue_ind = residue_ind == SIZE_MAX ? new_ind : std::max(residue_ind, new_ind);
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::incrResidueInd () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            residue_ind++;
        } else if ( fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z' ) {
            residue_ind--;
        } else { throw 30; }
    }

    bool ModifiedPeptide::FragmentGraph::isResidueEnd () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            return residue_ind == modified_peptide->residues.size();
        } else if ( fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z' ) {
            return residue_ind == SIZE_MAX;
        } else { throw 30; }
    }

    size_t ModifiedPeptide::FragmentGraph::getResidueDistance () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            return residue_ind;
        } else if ( fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z' ) {
            return modified_peptide->residues.size() - residue_ind - 1;
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::calculateFragment () {
        size_t mod_state = 0;
        if ( signature.count(residue_ind) ) {
            mod_state = signature.at(residue_ind);
        }

        float new_fragment_mz = modified_peptide->residues[residue_ind][mod_state];
        if (running_sum.size() != 0) {
            new_fragment_mz += running_sum.back();
        }
        running_sum.push_back(new_fragment_mz);

        running_sequence.push_back( modified_peptide->peptide[residue_ind] );
    }

    void ModifiedPeptide::FragmentGraph::updateLosses () {
        size_t mod_state = 0;
        if ( signature.count(residue_ind) ) {
            mod_state = signature.at(residue_ind);
        }

        if (modified_peptide->neutral_losses[residue_ind][mod_state] != 0.) {
            neutral_loss_stack.push_back( 
                modified_peptide->neutral_losses[residue_ind][mod_state] 
            );
            neutral_loss_iter.reset(neutral_loss_stack, 2);
        }
        running_loss_count.push_back( neutral_loss_stack.size() );
        neutral_loss_iter.reset(); // does nothing if already reset 
    }

    void ModifiedPeptide::FragmentGraph::resetIterator () {
	// First reset signature, then reset fragment
        modifiable.clear();
        signature.clear();

        n_mods_outstanding = modified_peptide->n_of_mod;
        for (resetResidueInd(); !isResidueEnd(); incrResidueInd()) {
            if ( isPositionModifiable() ) {
                modifiable.push_back(residue_ind);
                if ( n_mods_outstanding ) {
                    n_mods_outstanding--;
                    signature[modifiable.back()] = 1;
                } else {
                    signature[modifiable.back()] = 0;
                }
            }
        }

	resetFragment();
    }

    void ModifiedPeptide::FragmentGraph::incrSignature () {
        // Don't attempt to increment if signature end
        if (isSignatureEnd()) {throw 40;}
        
        // Otherwise, one mod is outstanding to start
        n_mods_outstanding = 1;
        size_t final_index = modified_peptide->residues.size();
        for (std::vector<size_t>::iterator it = --modifiable.end();
             it >= modifiable.begin() && n_mods_outstanding; 
             it--) {
            if ( signature.at(*it) > 0 ) {
                final_index = *it;
                signature.at(*it) = 0;

                size_t n_pos_until_end = modifiable.end() - it;
                if (  n_mods_outstanding == n_pos_until_end ) {
                    n_mods_outstanding++;
                    continue;
                }

                for (; n_mods_outstanding; n_mods_outstanding--) {
                    signature.at(*++it)++;
                }

            }
        }

        if (!isSignatureEnd()) { 
            setResidueInd(final_index);
            running_sum.resize(getResidueDistance());
            running_sequence.resize(getResidueDistance());

            running_loss_count.resize(getResidueDistance());
            if ( !running_loss_count.empty() ) {
                neutral_loss_stack.resize(running_loss_count.back());
            } else { neutral_loss_stack.clear(); }
            neutral_loss_iter.reset(neutral_loss_stack, 2);

            calculateFragment();
            updateLosses();
       }
    }

    bool ModifiedPeptide::FragmentGraph::isSignatureEnd () {
        return n_mods_outstanding > 0;
    }

    void ModifiedPeptide::FragmentGraph::resetFragment () {
        // Intialize running state trackers
        // Makes sure enough space is available and
        // clears away any old data.
        size_t peptide_size = modified_peptide->peptide.size();
        running_sum.reserve(peptide_size);
        running_sum.clear();

        running_sequence.reserve(peptide_size);
        running_sequence.clear();

        running_loss_count.reserve(peptide_size);
        running_loss_count.clear();
        neutral_loss_stack.clear();
        neutral_loss_iter.reset(neutral_loss_stack, 2);

        // Reset to first step
        resetResidueInd();
        calculateFragment();
        updateLosses();
    }

    void ModifiedPeptide::FragmentGraph::incrFragment () {
        // Don't increment if fragment end or signature end
        if (isSignatureEnd() || isFragmentEnd()) {throw 40;}

        // Only increment residue if more losses are not available
        if ( neutral_loss_iter.hasNext() ) {
            neutral_loss_iter.next();
        } else {
            incrResidueInd();
            if (!isFragmentEnd()) {
                calculateFragment();
                updateLosses();
            }
        }
    }

    bool ModifiedPeptide::FragmentGraph::isFragmentEnd () {
        bool last_fragment;
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            last_fragment = residue_ind == modified_peptide->residues.size() - 1;
        } else if ( fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z') {
            last_fragment = residue_ind == 0;
        } else { throw 30; }
        return last_fragment && !neutral_loss_iter.hasNext();
    }

    bool ModifiedPeptide::FragmentGraph::isLoss () {
        return neutral_loss_iter.getPos() > 0;
    }

    void ModifiedPeptide::FragmentGraph::setSignature (std::vector<size_t> new_signature) {
        if (new_signature.size() != modifiable.size()) {throw 50;}
        
        if (fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z') {
            std::reverse(new_signature.begin(), new_signature.end());
        }

        // Reset fragment state
        running_sum.clear();
        running_sequence.clear();
        running_loss_count.clear();
        neutral_loss_stack.clear();
        neutral_loss_iter.reset(neutral_loss_stack, 2);

        auto sig_it = new_signature.begin();
        auto mod_it = modifiable.begin();
        for (; sig_it != new_signature.end(); sig_it++, mod_it++){
            signature[*mod_it] = *sig_it;
        }    

        // Reset to first step
        resetResidueInd();
        calculateFragment();
        updateLosses();
    }

    std::vector<size_t> ModifiedPeptide::FragmentGraph::getSignature () {
        std::vector<size_t> signature_vector;
        for ( size_t it : modifiable ) {
            signature_vector.push_back(signature.at(it));
        }

        // For efficiency's sake, y fragment signatures are always backwards
        if (fragment_type == 'y' || fragment_type == 'z' || fragment_type == 'Z') {
            std::reverse(signature_vector.begin(), signature_vector.end());
        }

        return signature_vector;
    }

    float ModifiedPeptide::FragmentGraph::getFragmentMZ () {

        double fragment_mz = running_sum.back() - neutral_loss_iter.getSum();
        if (fragment_type == 'y') {
            fragment_mz += 18.010565; // Gain of H2O
        } else if (fragment_type == 'z') {
            fragment_mz += 18.010565; // Gain of H2O
            fragment_mz -= 17.026549; // Loss of NH3
        } else if (fragment_type == 'Z') {
            fragment_mz += 18.010565; // Gain of H2O
            fragment_mz -= 16.018724; // Loss of NH2
        } else if (fragment_type == 'c') {
            fragment_mz += 17.026549; // Gain of NH3
        }

        if (charge_state > 0) {
            fragment_mz = (fragment_mz + charge_state * 1.007825) / charge_state;
        }

        return fragment_mz;

    }

    size_t ModifiedPeptide::FragmentGraph::getFragmentSize () {
        return running_sequence.size();
    }

    std::string ModifiedPeptide::FragmentGraph::getFragmentSeq () {
        return running_sequence;
    }

    bool ModifiedPeptide::FragmentGraph::isPositionModifiable () {
        return modified_peptide->residues[residue_ind].size() > 1;
    }

    bool ModifiedPeptide::FragmentGraph::isPositionModified () {
        if ( signature.count(residue_ind) ) {
            return signature.at(residue_ind);
        } else { return false; }
    }

}
