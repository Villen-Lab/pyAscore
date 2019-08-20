#include <iostream>
#include <vector>
#include <string>
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

    void ModifiedPeptide::applyAuxMods (const size_t * aux_mod_pos = NULL,
                                        const float * aux_mod_mass = NULL,
                                        size_t n_aux_mods = 0) {

        for (size_t ind = 0; ind < n_aux_mods; ind++) {
            // A residue with an aux mod is not allowed to have more than 1 mod
            residues[aux_mod_pos[ind]].resize(1);
            residues[aux_mod_pos[ind]].front() += aux_mod_mass[ind];
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

    }

    void ModifiedPeptide::printInternals () {
        std::cout << std::endl << "Internal State:" << std::endl;
        for (size_t i = 0; i < residues.size(); i++) {
            std::cout << peptide[i];
            for (size_t j = 0; j < residues[i].size(); j++) {
                std::cout << "    " << residues[i][j];
            }
            std::cout << std::endl;
        }
    }

}
