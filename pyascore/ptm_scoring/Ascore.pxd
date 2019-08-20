cdef extern from "lib/Ascore.cpp":
    pass

cdef extern from "lib/Ascore.h" namespace "ptmscoring":
    cdef cppclass Ascore:
        Ascore() except +
        Ascore(float) except +
        float window_size
        void getWindowSize()
