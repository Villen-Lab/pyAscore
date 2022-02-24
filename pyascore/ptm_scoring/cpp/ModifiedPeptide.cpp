#include <iostream>
#include <cstdio>
#include <vector>
#include <tuple>
#include <string>
#include <algorithm>
#include <unordered_map>
#include <cmath>
#include <ctype.h>
#include <functional>
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
        for (char const& res : peptide) {
            residues.push_back({ StandardResidues.at(res) });

            // If character is in neutral loss group, add neutral loss
            neutral_losses.push_back({0.});
            if (nl_groups.count(res)) {
                neutral_losses.back().front() = nl_groups.at(res);
            } 

            // If character is in mod group, add a modified residue
            if (mod_group.find(res) != std::string::npos) {
                residues.back().push_back( StandardResidues.at(res) + mod_mass );
                neutral_losses.back().push_back(0.);
                if (nl_groups.count(std::tolower(res))) {
                    neutral_losses.back().back() = nl_groups.at(std::tolower(res));
                } 
            }

        }

    }

    void ModifiedPeptide::applyAuxMods () {

        for (size_t ind = 0; ind < aux_mod_pos.size(); ind++) {
            if (aux_mod_pos[ind] == 0) {
                transform(residues[aux_mod_pos[ind]].begin(), 
                          residues[aux_mod_pos[ind]].end(), 
                          residues[aux_mod_pos[ind]].begin(),
                          bind2nd(std::plus<float>(), aux_mod_mass[ind]));
            } else {
                size_t pep_ind = aux_mod_pos[ind] - 1;
                residues[pep_ind].front() += aux_mod_mass[ind];
                // A residue with a non n-terminus aux mod cannot have a variable mod
                residues[pep_ind].resize(1);
                
                char res = peptide.at(aux_mod_pos[ind] - 1);
                if (nl_groups.count(std::tolower(res))) {
                    neutral_losses[pep_ind].front() = nl_groups.at(std::tolower(res));
                    neutral_losses[pep_ind].resize(1);
                }
            }
        }
    }

    void ModifiedPeptide::initializeFragments () {
        fragment_scores.clear();
        fragments.clear();
        for (char t : fragment_types){
            // Only supporting charge 1 for now
            for (FragmentGraph graph = getFragmentGraph(t, 1);
                 !graph.isSignatureEnd();
                 graph.incrSignature()) {
                for(; !graph.isFragmentEnd();
                      graph.incrFragment()) {
                    fragments.push_back(graph.getFragmentMZ());
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

    void ModifiedPeptide::consumePeptide (std::string peptide, size_t n_of_mod,
                                          const unsigned int * aux_mod_pos,
                                          const float * aux_mod_mass,
                                          size_t n_aux_mods) {

        this->peptide = peptide;
        this->n_of_mod = n_of_mod;

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
            if (mz > ppm_low && mz < ppm_high) {
                if (!fragment_scores.count(*it) or std::get<1>(fragment_scores.at(*it)) > rank) {
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

    size_t ModifiedPeptide::getNumberModifiable () const {
        size_t n_modifiable = 0;
        for (const std::vector<float> & res_it : residues) {
            if (res_it.size() > 1) {n_modifiable += 1;} 
        }
        return n_modifiable;
    }

    size_t ModifiedPeptide::getPosOfNthModifiable (size_t n) const {
        for (size_t ind = 0; ind < residues.size(); ind++) {
            if (residues.at(ind).size() > 1 and n == 0) {
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
            signature = std::vector<size_t>(n_of_mod, 1);
        }

        std::string mod_peptide = "";

        // If an n terminal mod is present, pin it to the front
        size_t aux_ind = 0;
        if (!aux_mod_pos.empty() and aux_mod_pos.front() == 0) {
            mod_peptide += "n";

            char mod_buffer[10];
            std::sprintf(mod_buffer, "[%d]", (int) std::round(aux_mod_mass[0]));
            mod_peptide += mod_buffer;

            aux_ind++;
        }

        size_t sig_ind = 0;
        for (size_t pep_ind = 0; pep_ind < peptide.size(); pep_ind++) { 
            char aa = peptide[pep_ind];
            mod_peptide += aa;

            if (mod_group.find(aa) != std::string::npos 
                and sig_ind < signature.size() 
                and signature[sig_ind++] == 1) {
                // Add a modification mass
                char mod_buffer[10];
                std::sprintf(mod_buffer, "[%d]", (int) std::round(mod_mass));
                mod_peptide += mod_buffer;
            } else if (aux_ind < aux_mod_pos.size()
                       and (pep_ind + 1) == aux_mod_pos[aux_ind]) {
                // Add a mass specified in aux_masses
                char mod_buffer[10];
                std::sprintf(mod_buffer, "[%d]", (int) std::round(aux_mod_mass[aux_ind]));
                mod_peptide += mod_buffer;
                aux_ind++;
            }
        }
        return mod_peptide;
    }

    ModifiedPeptide::FragmentGraph ModifiedPeptide::getFragmentGraph (char fragment_type, size_t charge_state) const {
        return ModifiedPeptide::FragmentGraph(this, fragment_type, charge_state);
    }

    std::vector<std::vector<float>> ModifiedPeptide::getSiteDeterminingIons (const std::vector<size_t> & signature_1,
                                                                             const std::vector<size_t> & signature_2,
                                                                             char fragment_type, size_t charge_state) const {

        std::vector<std::tuple<float, size_t>> fragments;

        // Gather fragments from the first signature
        ModifiedPeptide::FragmentGraph graph_1 = getFragmentGraph(fragment_type, charge_state);
        graph_1.setSignature(signature_1);
        for(; !graph_1.isFragmentEnd();
               graph_1.incrFragment()) {
            fragments.push_back( {graph_1.getFragmentMZ(), 0} );
        }

        // Gather fragments from the second signature
        ModifiedPeptide::FragmentGraph graph_2 = getFragmentGraph(fragment_type, charge_state);
        graph_2.setSignature(signature_2);
        for(; !graph_2.isFragmentEnd();
               graph_2.incrFragment()) {
            fragments.push_back( {graph_2.getFragmentMZ(), 1} );
        }

        // Only store peaks that are different
        std::sort(fragments.begin(), fragments.end());
        std::vector<std::vector<float>> ion_lists(2);
        for (size_t frag_ind = 0; frag_ind < fragments.size(); frag_ind++) {
            bool keep = true;
            if ( frag_ind > 0 ) {
                keep = keep && std::abs(std::get<0>(fragments[frag_ind]) - std::get<0>(fragments[frag_ind - 1])) > mz_error;
            }
            if ( frag_ind + 1 < fragments.size() ) {
                keep = keep && std::abs(std::get<0>(fragments[frag_ind]) - std::get<0>(fragments[frag_ind + 1])) > mz_error;
            }
            if ( keep ) {
                ion_lists[ std::get<1>(fragments[frag_ind]) ].push_back( std::get<0>(fragments[frag_ind]) );
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
        } else if ( fragment_type == 'y'  || fragment_type == 'z' ) {
            residue_ind = modified_peptide->residues.size() - 1;
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::setResidueInd (size_t new_ind) {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            residue_ind = std::min(residue_ind, new_ind);
        } else if ( fragment_type == 'y' || fragment_type == 'z' ) {
            residue_ind = residue_ind == SIZE_MAX ? new_ind : std::max(residue_ind, new_ind);
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::incrResidueInd () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            residue_ind++;
        } else if ( fragment_type == 'y' || fragment_type == 'z' ) {
            residue_ind--;
        } else { throw 30; }
    }

    bool ModifiedPeptide::FragmentGraph::isResidueEnd () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            return residue_ind == modified_peptide->residues.size();
        } else if ( fragment_type == 'y' || fragment_type == 'z' ) {
            return residue_ind == SIZE_MAX;
        } else { throw 30; }
    }

    size_t ModifiedPeptide::FragmentGraph::getResidueDistance () {
        if ( fragment_type == 'b' || fragment_type == 'c' ) {
            return residue_ind;
        } else if ( fragment_type == 'y' || fragment_type == 'z' ) {
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

        // Reset to first step
        resetResidueInd();
        calculateFragment();
        updateLosses();
    }

    void ModifiedPeptide::FragmentGraph::incrSignature () {
        // Don't attempt to increment if signature end
        if (isSignatureEnd()) {throw 40;}
        
        // Otherwise, one mod is outstanding to start
        n_mods_outstanding = 1;
        size_t final_index = modified_peptide->residues.size();
        for (std::vector<size_t>::iterator it = --modifiable.end();
             it >= modifiable.begin() and n_mods_outstanding; 
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

    void ModifiedPeptide::FragmentGraph::incrFragment () {
        // Don't increment if fragment end or signature end
        if (isSignatureEnd() or isFragmentEnd()) {throw 40;}

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
        } else if ( fragment_type == 'y' || fragment_type == 'z' ) {
            last_fragment = residue_ind == 0;
        } else { throw 30; }
        return last_fragment && !neutral_loss_iter.hasNext();
    }

    bool ModifiedPeptide::FragmentGraph::isLoss () {
        return neutral_loss_iter.getPos() > 0;
    }

    void ModifiedPeptide::FragmentGraph::setSignature (std::vector<size_t> new_signature) {
        if (new_signature.size() != modifiable.size()) {throw 50;}
        
        if (fragment_type == 'y' || fragment_type == 'z') {
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
        if (fragment_type == 'y' || fragment_type == 'z') {
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
