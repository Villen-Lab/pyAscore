import unittest
import numpy as np
from pyascore import PyBinnedSpectra

class TestPyBinnedSpectra(unittest.TestCase):
    def test_init_dealoc(self):
        params = [dict(bin_size=100., n_top=10),
                  dict(bin_size=150., n_top=10)]

        for pars in params:
            spec = PyBinnedSpectra(**pars)
            self.assertEqual(spec.bin_size, pars["bin_size"])
            del spec

    def test_spectral_processing(self):
        masses = np.array([100., 300., 325., 350., 375., 400., 425., 450., 475., 500., 550., 1000.])
        intensities = np.array([50., 200., 100., 1000., 500., 100., 1000., 200., 300., 400., 500., 50.])

        true_n_peaks = iter([1, 6, 2, 1])
        true_rank_0 = iter([100., 425., 550., 1000.])

        spec = PyBinnedSpectra(bin_size=200., n_top=6)
        spec.consume_spectra(masses, intensities)
        
        self.assertEqual(spec.min_mz, 100.)
        self.assertEqual(spec.max_mz, 1000.)
        self.assertEqual(spec.n_bins, 5)

        while spec.bin < spec.n_bins:
            if spec.n_peaks > 0:
                self.assertEqual(spec.n_peaks, next(true_n_peaks))
                self.assertEqual(spec.mz, next(true_rank_0))

            spec.next_bin()
            spec.reset_rank();
            #self.assertEqual(spec.n_peaks, next(true_n_peaks))

            #spec.reset_rank()
            #self.assertEqual(spec.mz, true_rank_0[ind])

            #spec.next_rank()
            #self.assertEqual(spec.mz, true_rank_1[ind])
