from libcpp.string cimport string

cdef extern from "lib/ModifiedPeptide.cpp":
    pass

cdef extern from "lib/ModifiedPeptide.cpp" namespace "ptmscoring":
    cdef cppclass ModifiedPeptide:
        ModifiedPeptide(string, float);
        void consumePeptide(string, size_t)
        void printInternals();
