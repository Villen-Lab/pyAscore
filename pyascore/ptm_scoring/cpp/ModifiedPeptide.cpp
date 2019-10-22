#include <iostream>
#include <vector>
#include <tuple>
#include <string>
#include <algorithm>
#include <unordered_map>
#include "Types.h"
#include "ModifiedPeptide.h"

namespace ptmscoring {

    ModifiedPeptide::ModifiedPeptide (std::string mod_group, float mod_mass){
        this->mod_group = mod_group;
        this->mod_mass = mod_mass;
    }

    ModifiedPeptide::~ModifiedPeptide () {};

    void ModifiedPeptide::initializeResidues () {

        // Dealloc if already full
        residues.resize(0); 

        // Iterate through characters of peptide
        for (char res : peptide) {
            residues.push_back({ StandardResidues.at(res) });

            // If character is in mod group, add a modified residue
            if (mod_group.find(res) != std::string::npos) {
                residues.back().push_back( StandardResidues.at(res) + mod_mass );
            }

        }

    }

    void ModifiedPeptide::applyAuxMods (const size_t * aux_mod_pos,
                                        const float * aux_mod_mass,
                                        size_t n_aux_mods) {

        for (size_t ind = 0; ind < n_aux_mods; ind++) {
            residues[aux_mod_pos[ind]].front() += aux_mod_mass[ind];
            if ( residues[aux_mod_pos[ind]].size() > 1) {
                // A residue with an aux mod is not allowed to have more than 1 mod
                residues[aux_mod_pos[ind]].resize(1);
            }
        }

    }

    void ModifiedPeptide::initializeFragments () {
        std::vector<char> fragment_types = {'b', 'y'};
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

    void ModifiedPeptide::consumePeptide (std::string peptide, size_t n_of_mod,
                                          const size_t * aux_mod_pos,
                                          const float * aux_mod_mass,
                                          size_t n_aux_mods) {

        this->peptide = peptide;
        this->n_of_mod = n_of_mod;

        initializeResidues();
        applyAuxMods(aux_mod_pos, aux_mod_mass, n_aux_mods);

        fragment_scores.clear();
        fragments.clear();
        initializeFragments();

    }

    void ModifiedPeptide::consumePeak (float mz, size_t rank) {
        std::vector<float>::iterator it;
        it = lower_bound(fragments.begin(), fragments.end(), mz - .5);

        float ppm_low;
        float ppm_high;
        for(; it < fragments.end(); it++) {
            ppm_low = *it - .5;
            ppm_high = *it + .5;
            if (mz > ppm_low and mz < ppm_high) {
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
    
    std::string ModifiedPeptide::getModGroup() const {
        return mod_group;
    }

    float ModifiedPeptide::getModMass() const { 
        return mod_mass;
    } 

    std::string ModifiedPeptide::getPeptide() const {
        return peptide;
    }

    ModifiedPeptide::FragmentGraph ModifiedPeptide::getFragmentGraph (char fragment_type, size_t charge_state) const {
        return ModifiedPeptide::FragmentGraph(this, fragment_type, charge_state);
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
        if ( fragment_type == 'b' ) {
            residue_ind = 0;
        } else if ( fragment_type == 'y' ) {
            residue_ind = modified_peptide->residues.size() - 1;
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::setResidueInd (size_t new_ind) {
        if ( fragment_type == 'b' ) {
            residue_ind = std::min(residue_ind, new_ind);
        } else if ( fragment_type == 'y' ) {
            residue_ind = residue_ind == SIZE_MAX ? new_ind : std::max(residue_ind, new_ind);
        } else { throw 30; }
    }

    void ModifiedPeptide::FragmentGraph::incrResidueInd () {
        if ( fragment_type == 'b' ) {
            residue_ind++;
        } else if ( fragment_type == 'y' ) {
            residue_ind--;
        } else { throw 30; }
    }

    bool ModifiedPeptide::FragmentGraph::isResidueEnd () {
        if ( fragment_type == 'b' ) {
            return residue_ind == modified_peptide->residues.size();
        } else if ( fragment_type == 'y' ) {
            return residue_ind == SIZE_MAX;
        } else { throw 30; }
    }

    size_t ModifiedPeptide::FragmentGraph::getResidueDistance () {
        if ( fragment_type == 'b' ) {
            return residue_ind;
        } else if ( fragment_type == 'y' ) {
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

    void ModifiedPeptide::FragmentGraph::resetIterator () {
        // Intialize running state trackers
        // Makes sure enough space is available and
        // clears away any old data.
        size_t peptide_size = modified_peptide->peptide.size();
        running_sum.reserve(peptide_size);
        running_sum.clear();

        running_sequence.reserve(peptide_size);
        running_sequence.clear();

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
            calculateFragment();
       }
    }

    bool ModifiedPeptide::FragmentGraph::isSignatureEnd () {
        return n_mods_outstanding > 0;
    }

    void ModifiedPeptide::FragmentGraph::incrFragment () {
        // Don't increment if fragment end or signature end
        if (isSignatureEnd() or isFragmentEnd()) {throw 40;}

        // Otherwise, increment the residue ind
        incrResidueInd();
        if (!isFragmentEnd()) {calculateFragment();};
    }

    bool ModifiedPeptide::FragmentGraph::isFragmentEnd () {
        return isResidueEnd();
    }

    std::vector<size_t> ModifiedPeptide::FragmentGraph::getSignature () {
        std::vector<size_t> signature_vector;
        for ( size_t it : modifiable ) {
            signature_vector.push_back(signature.at(it));
        }

        // For efficiency's sake, y fragment signatures are always backwards
        if (fragment_type == 'y') {
            std::reverse(signature_vector.begin(), signature_vector.end());
        }

        return signature_vector;
    }

    float ModifiedPeptide::FragmentGraph::getFragmentMZ () {

        double fragment_mz = running_sum.back();
        if (fragment_type == 'y') {
            fragment_mz += 18.01528;
        }

        if (charge_state > 0) {
            fragment_mz = (fragment_mz + charge_state) / charge_state;
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
