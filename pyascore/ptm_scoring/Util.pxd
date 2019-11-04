cdef extern from "cpp/Util.cpp":
    pass

cdef extern from "cpp/Util.cpp" namespace "ptmscoring":
    cdef cppclass LogMath:
        float log_sum(float, float);
        float log_bin_coef(size_t, size_t);

    cdef cppclass BinomialDist:
        BinomialDist(float);
        float log_pmf(size_t, size_t);
        float log_pvalue(size_t, size_t);
        float log10_pvalue(size_t, size_t);
