#include <iostream>
#include "Ascore.h"
#include "Spectra.h"
#include "ModifiedPeptide.h"

namespace ptmscoring {

    Ascore::Ascore () {}

    Ascore::~Ascore () {}

    void Ascore::score (const BinnedSpectra & spectra, const ModifiedPeptide & peptide) {
        std::cout << std::endl << "Scoring not implemented yet" << std::endl;
    }
}
