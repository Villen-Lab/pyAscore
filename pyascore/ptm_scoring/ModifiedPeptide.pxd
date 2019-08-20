from libcpp.string cimport string

cdef extern from "cpp/ModifiedPeptide.cpp":
    pass

cdef extern from "cpp/ModifiedPeptide.cpp" namespace "ptmscoring":
    cdef cppclass ModifiedPeptide:
        ModifiedPeptide(string, float);
        void consumePeptide(string, size_t)
