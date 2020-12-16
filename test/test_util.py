import unittest
import numpy as np
from scipy.special import logsumexp
from scipy.special import binom as binom_coef
from scipy.stats import binom as binom_dist
from pyascore import PyLogMath, PyBinomialDist, PyPowerSetSum

class TestLogMath(unittest.TestCase):
    lmath = PyLogMath()

    def test_log_sum(self):
        summands = [(-np.inf, 0),
                    (0, -np.inf)]

        n_to_add = 100
        np.random.seed(2345)
        for a, b in zip(np.random.randn(n_to_add), np.random.randn(n_to_add)):
            summands.append((a, b))

        for s in summands:
            self.assertTrue(np.isclose(self.lmath.log_sum(*s), 
                                       logsumexp(s),
                                       rtol=0, atol=1e-6))

    def test_log_bin_coef(self):
        MAX_N = 50

        for n in range(1, MAX_N + 1):
            for k in range(1, n + 1):
                #print(k, n, self.lmath.log_bin_coef(k, n), np.log(binom_coef(n, k)))
                self.assertTrue(np.isclose(self.lmath.log_bin_coef(k, n),
                                           np.log(binom_coef(n, k)),
                                           rtol=0, atol=5e-5))

class TestBinomDist(unittest.TestCase):
    def test_log_pdf(self):
        MAX_N = 50

        probs = [.1, .25, .5, .75, .9]
        dists = [PyBinomialDist(p) for p in probs]

        for p, d in zip(probs, dists):
            for n in range(MAX_N):
                for k in range(1, n + 1):
                    #print(k, n, d.log_pmf(k, n), binom_dist.logpmf(k, n, p))
                    self.assertTrue(np.isclose(d.log_pmf(k, n),
                                               binom_dist.logpmf(k, n, p),
                                               rtol=0, atol=5e-5))
    def test_log_pvalue(self):
        MAX_N = 50

        probs = [.1, .25, .5, .75, .9]
        dists = [PyBinomialDist(p) for p in probs]

        for p, d in zip(probs, dists):
            for n in range(MAX_N):
                for k in range(1, n + 1):
                    scipy_pvalue = logsumexp([binom_dist.logpmf(k, n, p),
                                              binom_dist.logsf(k, n, p)])
                    #print(k, n, d.log_pvalue(k, n), scipy_pvalue)
                    self.assertTrue(np.isclose(d.log_pvalue(k, n),
                                               scipy_pvalue,
                                               rtol=0, atol=5e-5))
                    self.assertTrue(np.isclose(d.log10_pvalue(k, n),
                                               np.log10(np.exp(scipy_pvalue)),
                                               rtol=0, atol=5e-5))

class TestPowerSetSum(unittest.TestCase):
    def test_empty_set(self):
        pss = PyPowerSetSum()

        self.assertFalse( pss.has_next() )
        self.assertTrue( pss.get_sum() == 0.)

    def test_simple_sums(self):
        target = np.array([1., 2., 3.], dtype=np.float32)
        pss = PyPowerSetSum(target, 2)

        self.assertTrue( pss.has_next() )
        self.assertTrue( pss.get_sum() == 0.)

        true_sums = iter([1., 2., 3., 4., 5.])
        while pss.has_next():
            pss.next()
            self.assertTrue( pss.get_sum() == next(true_sums) )

        new_target = np.array([4., 5., 6.], dtype=np.float32)
        pss.reset(new_target, 2)

        self.assertTrue( pss.has_next() )
        self.assertTrue( pss.get_sum() == 0.)

        true_sums = iter([4., 5., 6., 9., 10, 11.])
        while pss.has_next():
            pss.next()
            self.assertTrue( pss.get_sum() == next(true_sums) )
