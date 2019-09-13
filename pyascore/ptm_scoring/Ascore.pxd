from ModifiedPeptide cimport ModifiedPeptide
from Spectra cimport BinnedSpectra

cdef extern from "cpp/Ascore.cpp":
    pass

cdef extern from "cpp/Ascore.h" namespace "ptmscoring":
    cdef cppclass Ascore:
        Ascore() except +
        void score(const BinnedSpectra &, const ModifiedPeptide &)
