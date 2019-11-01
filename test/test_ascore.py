import unittest
import time
import os
import pickle
import numpy as np
from pyascore import PyAscore

class TestPyAscore(unittest.TestCase):
    def test_single_spectrum_score(self):
        ascore = PyAscore(min_mz=500., max_mz=1500., bin_size=100., n_top=10,
                          mod_group="STY", mod_mass=79.966331)
        def test_match_pairs(match_file, spec_file):
            with open(os.path.join("test", "match_spectra_pairs", match_file), "rb") as src:
                match_list = pickle.load(src)

            with open(os.path.join("test", "match_spectra_pairs", spec_file), "rb") as src:
                spectra_list = pickle.load(src)

            n = 0
            running_time = 0
            for match, spectra in zip(match_list, spectra_list):
                #print(match)
                running_time -= time.time()
                ascore = PyAscore(min_mz=500., max_mz=1500., bin_size=100., n_top=10,
                                  mod_group="STY", mod_mass=79.966331)
                ascore.score(spectra["mz_values"], spectra["intensity_values"], match["peptide"], len(match["mod_positions"]))
                running_time += time.time()
                #print(ascore.best_sequence, ascore.best_score)
                n += 1

            print("Average for ({}, {}): {} sec/match".format( match_file, spec_file, (running_time/n) ))
            print()

        test_match_pairs("velos_matches_3_mods.pkl", "velos_spectra_3_mods.pkl")
        test_match_pairs("dump_match.pkl", "dump_spectra.pkl")
