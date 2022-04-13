import numpy as np
cimport numpy as np
from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "cpp/ModifiedPeptide.cpp":
    pass

cdef extern from "cpp/ModifiedPeptide.cpp" namespace "ptmscoring":
    cdef cppclass ModifiedPeptide:
        ModifiedPeptide(string, float, float, string)
        void addNeutralLoss(string, float)
        void consumePeptide(string, size_t, size_t)
        void consumePeptide(string, size_t, size_t,
                            const unsigned int *, 
                            const float *,
                            size_t)
        string getPeptide(vector[size_t])
        void consumePeak(float, size_t)
        size_t getNumberOfMods()

        cppclass FragmentGraph:
            FragmentGraph(const ModifiedPeptide *, char, size_t);

            char getFragmentType();
            size_t getChargeState();

            void resetIterator();
            void incrSignature();
            bint isSignatureEnd();
            void resetFragment();
            void incrFragment();
            bint isFragmentEnd();

            void setSignature(vector[size_t]);
            vector[size_t] getSignature();
            float getFragmentMZ();
            size_t getFragmentSize();
            string getFragmentSeq();

        FragmentGraph getFragmentGraph(char, size_t);
        vector[vector[float]] getSiteDeterminingIons(
            const vector[size_t] &,
            const vector[size_t] &,
            char, size_t
        )
