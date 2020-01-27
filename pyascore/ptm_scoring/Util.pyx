# distutils: language = c++
import cython

from Util cimport LogMath, BinomialDist

cdef class PyLogMath:
    """
    The PyLogMath object provides access to efficient log based mathematical operations.
    """
    cdef LogMath * log_math_ptr

    def __cinit__(self):
        self.log_math_ptr = new LogMath()

    def __dealloc__(self):
        del self.log_math_ptr

    def log_sum(self, float a, float b):
        """Exponentiate arguments, sum, and log results.

        Computes:
            log(exp(a) + exp(b))

        Parameters
        ----------
        a : float
            Argument 1.
        b : float
            Argument 2.
        """
        return self.log_math_ptr[0].log_sum(a, b)

    def log_bin_coef(self, size_t k, size_t n):
        """Calculate the binomial coefficent in log space.

        Computes:
            log( n choose k )

        Parameters
        ----------
        k : size_t
            Subset size (number of successes).
        n : size_t
            Superset size (number of trials).
        """
        return self.log_math_ptr[0].log_bin_coef(k, n)

cdef class PyBinomialDist:
    """
    The PyLogMath object provides a cached version of the Binomial distribution.

    The caching of previously calculated distributional statistics leads to substantial
    speedups in scoring. This also means that a new object must be created anytime a new
    probability is desired. However, since most analyses use a single probability for a
    whole file, this is no problem.

    Parameters
    ----------
    prob : float
        Probability of success.
    """
    cdef BinomialDist * binomial_dist_ptr

    def __cinit__(self, float prob):
        self.binomial_dist_ptr = new BinomialDist(prob)

    def __dealloc__(self):
        del self.binomial_dist_ptr

    def log_pmf(self, size_t successes, size_t trials):
        """Log Probability Mass Function for the Binomial Distribution

        Parameters
        ----------
        successes : size_t
        trials : size_t
        """
        return self.binomial_dist_ptr[0].log_pmf(successes, trials)

    def log_pvalue(self, size_t successes, size_t trials):
        """The log of the probability of observing at least as many successes

        Parameters
        ----------
        successes : size_t
        trials : size_t
        """
        return self.binomial_dist_ptr[0].log_pvalue(successes, trials)

    def log10_pvalue(self, size_t successes, size_t trials):
        """The log10 of the probability of observing at least as many successes

        Parameters
        ----------
        successes : size_t
        trials : size_t
        """
        return self.binomial_dist_ptr[0].log10_pvalue(successes, trials)

