# distutils: language = c++
import cython

from Util cimport LogMath, BinomialDist

cdef class PyLogMath:
    cdef LogMath * log_math_ptr

    def __cinit__(self):
        self.log_math_ptr = new LogMath()

    def __dealloc__(self):
        del self.log_math_ptr

    def log_sum(self, float a, float b):
        return self.log_math_ptr[0].log_sum(a, b)

    def log_bin_coef(self, size_t k, size_t n):
        return self.log_math_ptr[0].log_bin_coef(k, n)

cdef class PyBinomialDist:
    cdef BinomialDist * binomial_dist_ptr

    def __cinit__(self, float prob):
        self.binomial_dist_ptr = new BinomialDist(prob)

    def __dealloc__(self):
        del self.binomial_dist_ptr

    def log_pmf(self, size_t successes, size_t trials):
        return self.binomial_dist_ptr[0].log_pmf(successes, trials)

    def log_pvalue(self, size_t successes, size_t trials):
        return self.binomial_dist_ptr[0].log_pvalue(successes, trials)
