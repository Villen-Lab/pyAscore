from libcpp.vector cimport vector
from libcpp.string cimport string

cdef extern from "cpp/ModifiedPeptide.cpp":
    pass

cdef extern from "cpp/ModifiedPeptide.cpp" namespace "ptmscoring":
    cdef cppclass ModifiedPeptide:
        ModifiedPeptide(string, float);
        void consumePeptide(string, size_t)

        void resetIterator(char);
        size_t incrSignature();
        vector[size_t] getSignature();
        size_t incrFragment();
        char getFragmentType();
        float getFragmentMZ(size_t);
        size_t getFragmentSize();
        string getFragmentSeq();
