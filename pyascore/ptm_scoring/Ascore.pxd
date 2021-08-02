from ModifiedPeptide cimport ModifiedPeptide
from Spectra cimport BinnedSpectra
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "cpp/Ascore.cpp":
    pass

cdef extern from "cpp/Ascore.h" namespace "ptmscoring":
    cdef struct ScoreContainer:
        vector[size_t] signature;
        vector[size_t] counts;
        vector[float] scores;
        float weighted_score;
        size_t total_fragments;

    cdef cppclass Ascore:
        Ascore() except +;
        void score(const BinnedSpectra &, const ModifiedPeptide &);
        float calculateAmbiguity(const ScoreContainer&, const ScoreContainer&);
        string getBestSequence();
        vector[string] getAllSequences();
        float getBestScore();
        vector[ScoreContainer] getAllPepScores();
        vector[float] getAscores();
        vector[size_t] getAlternativeSites(size_t);
