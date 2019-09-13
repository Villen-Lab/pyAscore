from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "cpp/ModifiedPeptide.cpp":
    pass

cdef extern from "cpp/ModifiedPeptide.cpp" namespace "ptmscoring":
    cdef cppclass ModifiedPeptide:
        ModifiedPeptide(string, float);
        void consumePeptide(string, size_t)
        void consumePeak(float, size_t);

        cppclass FragmentGraph:
            FragmentGraph(const ModifiedPeptide *, char, size_t);

            char getFragmentType();
            size_t getChargeState();

            void resetIterator();
            void incrSignature();
            bint isSignatureEnd();
            void incrFragment();
            bint isFragmentEnd();

            vector[size_t] getSignature();
            float getFragmentMZ();
            size_t getFragmentSize();
            string getFragmentSeq();

        FragmentGraph getFragmentGraph(char, size_t);
