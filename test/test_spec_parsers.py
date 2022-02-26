import unittest
import numpy as np
from pyascore import spec_parsers

SCAN_NUMBERS = [14760, 18330, 20462, 21996, 26219, 26962, 27845, 31328, 32257, 35669]
PRECURSOR_MZ = [846.306451825194, 871.696163579367, 858.378601074219, 1116.095703125,
                858.408142089844, 1427.79736328125, 827.992004394531, 1078.430162374319, 
                1023.712707519531, 885.028560474252]

class TestSpectraParser(unittest.TestCase):

    def test_mzml_reader(self):
        target_file="test/example_inputs/spectra/test_spectra.mzML"
        parser = spec_parsers.SpectraParser(target_file, "mzML")
        spectra_list = parser.to_list()

        self.assertEqual(len(spectra_list), len(SCAN_NUMBERS))
        for ind in range(len(spectra_list)):
            self.assertEqual(spectra_list[ind]["scan"], SCAN_NUMBERS[ind])
            self.assertEqual(spectra_list[ind]["ms_level"], 2)
            self.assertEqual(spectra_list[ind]["precursor_mz"], PRECURSOR_MZ[ind])
            self.assertEqual(spectra_list[ind]["precursor_charge"], 3)
            self.assertTrue(spectra_list[ind]["mz_values"].shape[0] > 0)
            self.assertTrue(spectra_list[ind]["intensity_values"].shape[0] > 0)

    def test_mzxml_reader(self):
        target_file="test/example_inputs/spectra/test_spectra.mzXML"
        parser = spec_parsers.SpectraParser(target_file, "mzXML")
        spectra_list = parser.to_list()

        self.assertEqual(len(spectra_list), len(SCAN_NUMBERS))
        for ind in range(len(spectra_list)):
            self.assertEqual(spectra_list[ind]["scan"], SCAN_NUMBERS[ind])
            self.assertEqual(spectra_list[ind]["ms_level"], 2)
            self.assertEqual(spectra_list[ind]["precursor_mz"], PRECURSOR_MZ[ind])
            self.assertEqual(spectra_list[ind]["precursor_charge"], 3)
            self.assertTrue(spectra_list[ind]["mz_values"].shape[0] > 0)
            self.assertTrue(spectra_list[ind]["intensity_values"].shape[0] > 0)
