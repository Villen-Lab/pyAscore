import unittest
import numpy as np
from pyascore import PyBinnedSpectra

class TestPyBinnedSpectra(unittest.TestCase):
    def test_init_dealoc(self):
        params = [dict(min_mz=500., max_mz=1500., bin_size=100., n_top=10),
                  dict(min_mz=500., max_mz=1500., bin_size=150., n_top=10)]

        n_bins_true = [10, 7]

        for pars, bins in zip(params, n_bins_true):
            spec = PyBinnedSpectra(**pars)
            self.assertEqual(spec.min_mz, pars["min_mz"])
            self.assertEqual(spec.max_mz, pars["max_mz"])
            self.assertEqual(spec.bin_size, pars["bin_size"])
            self.assertEqual(spec.n_bins, bins)
            del spec

    def test_spectral_processing(self):
        masses = np.array([100., 300., 325., 350., 375., 400., 425., 450., 475., 500., 550., 1000.])
        intensities = np.array([50., 200., 100., 1000., 500., 100., 1000., 200., 300., 400., 500., 50.])

        true_n_peaks = [3, 3, 2]
        true_rank_0 = [350., 425., 550.]
        true_rank_1 = [375., 475., 500.]

        spec = PyBinnedSpectra(min_mz=300., max_mz=600., bin_size=100., n_top=3)
        spec.consume_spectra(masses, intensities)
        for ind in range(spec.n_bins):
            spec.bin = ind
            self.assertEqual(spec.n_peaks, true_n_peaks[ind])

            spec.reset_rank()
            self.assertEqual(spec.mz, true_rank_0[ind])

            spec.next_rank()
            self.assertEqual(spec.mz, true_rank_1[ind])
