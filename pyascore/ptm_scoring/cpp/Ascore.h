#ifndef ASCORE_H
#define ASCORE_H

#include "Spectra.h"
#include "ModifiedPeptide.h"

namespace ptmscoring {

    class Ascore {
        public:
            Ascore();
            ~Ascore();

            void score(const BinnedSpectra &, const ModifiedPeptide &);
    };

}

#endif
