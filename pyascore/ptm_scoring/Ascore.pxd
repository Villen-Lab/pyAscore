cdef extern from "cpp/Ascore.cpp":
    pass

cdef extern from "cpp/Ascore.h" namespace "ptmscoring":
    cdef cppclass Ascore:
        Ascore() except +
        Ascore(float) except +
        float window_size
        void getWindowSize()
