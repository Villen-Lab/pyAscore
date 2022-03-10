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
        intensities = np.array([50., 200., 100., 1000., 500., 100., 1200., 200., 300., 400., 500., 50.])

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
            spec.reset_rank()

    def test_full_spectra_parse(self):
        # Build random spectra
        n_top = 10
        bin_size = 100.
        n_peaks = 500
        np.random.seed(2345)
        masses = np.random.uniform(500., 2000., n_peaks)
        intensities = 100. * np.random.randn(n_peaks) + 300.

        mz_low = np.min(masses)
        mz_high = np.max(masses)
        n_bins = int(np.ceil( ( mz_high - mz_low ) / bin_size ))
        
        spec = PyBinnedSpectra(bin_size=100., n_top=n_top)
        spec.consume_spectra(masses, intensities)

        for ind in range(n_bins):
            select = np.logical_and(masses >= (mz_low + ind * bin_size),
                                    masses < (mz_low + (ind + 1) * bin_size))
            binned_masses = masses[select]
            binned_intensities = intensities[select]

            sorted_ind = np.argsort(binned_intensities)[::-1]
            binned_masses = binned_masses[sorted_ind]
            binned_intensities = binned_intensities[sorted_ind]
            for rank in range(min(n_top, len(binned_masses))):
                self.assertEqual(spec.mz, binned_masses[rank])
                self.assertEqual(spec.intensity, binned_intensities[rank])
                spec.next_rank()
            spec.next_bin()
            spec.reset_rank()
