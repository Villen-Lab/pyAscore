#ifndef MODIFIEDPEPTIDE_H
#define MODIFIEDPEPTIDE_H

#include <vector>
#include <string>
#include <tuple>
#include <unordered_map>
#include "Types.h"
#include "Util.h"

namespace ptmscoring {

    class ModifiedPeptide {
        std::string mod_group;
        float mod_mass;
        float mz_error;
        std::string fragment_types;

        std::unordered_map<char, float> nl_groups;

        std::string peptide;
        size_t n_of_mod;
	size_t max_fragment_charge;
        std::vector<unsigned int> aux_mod_pos;
        std::vector<float> aux_mod_mass;

        std::vector<std::vector<float>> residues;
        std::vector<std::vector<float>> neutral_losses;
        std::vector<float> fragments;
        std::unordered_map<float, std::tuple<float, size_t>> fragment_scores;

        void initializeResidues();
        void applyAuxMods();
        void initializeFragments();
        public:
            ModifiedPeptide(std::string, float, float, std::string);
            ~ModifiedPeptide();

            void addNeutralLoss(std::string, float);

            void consumePeptide(std::string, size_t, 
			        size_t = 1,
                                const unsigned int * = NULL, 
                                const float * = NULL, 
                                size_t = 0);
            void consumePeak(float, size_t);
            bool hasMatch(float) const;
            std::tuple<float, size_t> getMatch(float) const;

            std::string getModGroup() const;
            float getModMass() const;
            float getMZError() const;
            std::string getFragmentTypes() const;
            size_t getNumberOfMods() const;
	    size_t getMaxFragmentCharge() const;
            size_t getNumberModifiable() const;
            size_t getPosOfNthModifiable(size_t) const;
            std::string getBasePeptide() const;
            std::string getPeptide(std::vector<size_t> = {}) const;

            class FragmentGraph;
            FragmentGraph getFragmentGraph(char, size_t) const;
            std::vector<std::vector<float>> getSiteDeterminingIons(
                const std::vector<size_t> &, 
                const std::vector<size_t> &,
                char, size_t
            ) const;
    };

    class ModifiedPeptide::FragmentGraph {
        const ModifiedPeptide * modified_peptide;
        size_t residue_ind;
        void resetResidueInd();
        void setResidueInd(size_t);
        void incrResidueInd();
        bool isResidueEnd();
        size_t getResidueDistance();

        char fragment_type;
        size_t charge_state;

        size_t n_mods_outstanding;
        std::vector<size_t> modifiable;
        std::unordered_map<size_t, size_t> signature;

        std::vector<float> running_sum;
        std::string running_sequence;
        std::vector<size_t> running_loss_count;
        std::vector<float> neutral_loss_stack;
        PowerSetSum neutral_loss_iter;
        void calculateFragment();
        void updateLosses();
        public:
            FragmentGraph(const ModifiedPeptide *, char, size_t);
            ~FragmentGraph();

            char getFragmentType();
            size_t getChargeState();

            void resetIterator();
            void incrSignature();
            bool isSignatureEnd();
	    void resetFragment();
            void incrFragment();
            bool isFragmentEnd();
            bool isLoss();

            void setSignature(std::vector<size_t>);
            std::vector<size_t> getSignature();
            float getFragmentMZ();
            size_t getFragmentSize();
            std::string getFragmentSeq();

            bool isPositionModifiable();
            bool isPositionModified();
    };


}

#endif
