#include <iostream>
#include <vector>
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

    void ModifiedPeptide::consumePeptide (std::string peptide, size_t n_of_mod,
                                          const size_t * aux_mod_pos,
                                          const float * aux_mod_mass,
                                          size_t n_aux_mods) {

        this->peptide = peptide;
        this->n_of_mod = n_of_mod;
        initializeResidues();
        applyAuxMods(aux_mod_pos, aux_mod_mass, n_aux_mods);
        resetIterator('b');

    }
    
    void ModifiedPeptide::resetIterator (char fragment_type){
        // Fragment type can only be changed with reinitialization
        this->fragment_type = fragment_type;

        // Intialize running state trackers
        // Makes sure enough space is available and
        // clears away any old data.
        running_sum.reserve(peptide.size() + 1);
        running_sum.clear();
        running_sum.push_back(0.);

        running_sequence.reserve(peptide.size());
        running_sequence.clear();

        modifiable.clear();
        signature.clear();
        size_t n_to_modify = n_of_mod;
        for (size_t ind = 0; ind < residues.size(); ind++) {
            if ( residues[ind].size() > 1 ) {
                modifiable.push_back(ind);
                if ( n_to_modify ) {
                    n_to_modify--;
                    signature[modifiable.back()] = 1;
                } else {
                    signature[modifiable.back()] = 0;
                }
            }
        }

        // Set to first fragment
        fragment_ind = -1;
        incrFragment();
    }

    size_t ModifiedPeptide::incrSignature () {
        size_t final_index = residues.size();
        size_t n_to_incr = 1;
        std::vector<size_t>::iterator it = modifiable.end();
        while ( it != modifiable.begin() and n_to_incr ) {
            it--;
            if ( signature.at(*it) > 0 ) {
                final_index = *it;
                signature.at(*it) = 0;
                
                size_t n_pos_until_end = modifiable.end() - it;
                if (  n_to_incr == n_pos_until_end ) {
                    n_to_incr++;
                    continue;
                }

                while ( n_to_incr ) {
                    it++;
                    signature.at(*it)++;
                    n_to_incr--;
                }

            }
        }

        if (n_to_incr) {
            return 0;
        } else {
            running_sum.resize(final_index + 1);
            running_sequence.resize(final_index);
            fragment_ind = final_index - 1;
            incrFragment();
            return 1;
       }
    }

    std::vector<size_t> ModifiedPeptide::getSignature() {

        std::vector<size_t> signature_vector;
        for ( size_t it : modifiable ) {
            signature_vector.push_back(signature.at(it));
        }

        return signature_vector;

    }

    void ModifiedPeptide::calculateFragment(size_t mod_state) {
        running_sequence.push_back( peptide[fragment_ind] );
        running_sum.push_back( running_sum.back() + residues[fragment_ind][mod_state] );
    }

    size_t ModifiedPeptide::incrFragment () {

        fragment_ind++;
        if ( fragment_ind >= peptide.size() ) {
            return 0;
        }

        if ( signature.count(fragment_ind) ) {
            calculateFragment(signature.at(fragment_ind));
        } else {
            calculateFragment(0);
        }

        return 1;

    }

    char ModifiedPeptide::getFragmentType () {
        return fragment_type;
    }

    float ModifiedPeptide::getFragmentMZ (size_t charge) {

        double fragment_mz = running_sum.back();
        if (fragment_type == 'y') {
            fragment_mz += 18.01528;
        }
        return (fragment_mz + charge) / charge;

    }

    size_t ModifiedPeptide::getFragmentSize () {
        return running_sequence.size();
    }

    std::string ModifiedPeptide::getFragmentSeq () {
        return running_sequence;
    }

}
