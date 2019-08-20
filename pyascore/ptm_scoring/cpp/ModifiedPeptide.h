#ifndef MODIFIEDPEPTIDE_H
#define MODIFIEDPEPTIDE_H

#include <vector>
#include <string>
#include <unordered_map>
#include "Types.h"

namespace ptmscoring {

    class ModifiedPeptide {
        std::string mod_group;
        float mod_mass;

        std::string peptide;
        size_t n_of_mod;

        std::vector<std::vector<float>> residues;
        std::vector<float> fragments;
        std::unordered_map<float, std::vector<size_t>> fragment_scores;

        void initializeResidues();
        void applyAuxMods(const size_t* = NULL, const float* = NULL, size_t = 0);
        void initializeFragments();
        public:
            ModifiedPeptide(std::string, float);
            ~ModifiedPeptide();

            void consumePeptide(std::string, size_t, const size_t * = NULL, const float * = NULL, size_t = 0);
            void consumePeak(Peak, size_t);

            std::string getModGroups(size_t) const;
            float getModMass(size_t) const;
    };

}

#endif
