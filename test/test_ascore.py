import unittest
import time
import os
import re
import pickle
import numpy as np
from pyascore import PyAscore

class TestPyAscore(unittest.TestCase):
    def test_single_spectrum_score(self):
        ascore = PyAscore(bin_size=100., n_top=10,
                          mod_group="STY", mod_mass=79.966331)
        def test_match_pairs(match_file, spec_file):
            with open(os.path.join("test", "match_spectra_pairs", match_file), "rb") as src:
                match_list = pickle.load(src)

            with open(os.path.join("test", "match_spectra_pairs", spec_file), "rb") as src:
                spectra_list = pickle.load(src)

            n = 0
            running_time = 0
            for match, spectra in zip(match_list, spectra_list):
                running_time -= time.time()
                ascore = PyAscore(bin_size=100., n_top=10,
                                  mod_group="STY", mod_mass=79.966331)
                ascore.score(spectra["mz_values"], 
                             spectra["intensity_values"], 
                             match["peptide"], 
                             len(match["mod_positions"]),
                             match["charge_state"] - 1)
                running_time += time.time()
                n += 1

            print("Average for ({}, {}): {} sec/match".format( 
                match_file, spec_file, (running_time/n) )
            )

            n = 0
            running_time = 0
            for match, spectra in zip(match_list, spectra_list):
                running_time -= time.time()
                ascore = PyAscore(bin_size=100., n_top=10,
                                  mod_group="STY", mod_mass=79.966331)
                ascore.add_neutral_loss("ST", 18.01528)
                ascore.score(spectra["mz_values"],
                             spectra["intensity_values"],
                             match["peptide"],
                             len(match["mod_positions"]),
                             match["charge_state"] - 1)
                running_time += time.time()
                n += 1

            print("Average with neutral loss for ({}, {}): {} sec/match".format(
                match_file, spec_file, (running_time/n) )
            )
            print()

        test_match_pairs("velos_matches_1_mods.pkl", "velos_spectra_1_mods.pkl")
        test_match_pairs("velos_matches_2_mods.pkl", "velos_spectra_2_mods.pkl")
        test_match_pairs("velos_matches_3_mods.pkl", "velos_spectra_3_mods.pkl")
        test_match_pairs("velos_matches_aux.pkl", "velos_spectra_aux.pkl")

    def test_alternative_site_consistency(self):
        ascore = PyAscore(bin_size=100., n_top=10,
                          mod_group="STY", mod_mass=79.966331)
        def test_consistency(match_file, spec_file):
            with open(os.path.join("test", "match_spectra_pairs", match_file), "rb") as src:
                match_list = pickle.load(src)


            with open(os.path.join("test", "match_spectra_pairs", spec_file), "rb") as src:
                spectra_list = pickle.load(src)

            for match, spectra in zip(match_list, spectra_list):
                ascore = PyAscore(bin_size=100., n_top=10,
                                  mod_group="STY", mod_mass=79.966331)
                ascore.score(spectra["mz_values"],
                             spectra["intensity_values"],
                             match["peptide"],
                             len(match["mod_positions"]))

                # Alternative sites should be free of duplicates
                for alt in ascore.alt_sites:
                    n_alt_sites = alt.shape[0]
                    n_deduplicated = np.unique(alt).shape[0]
                    self.assertEqual(n_alt_sites, n_deduplicated)

                # No modified sites should show up in alternative sites
                site_iter = re.finditer("[A-Z][^A-Z]*", ascore.best_sequence)
                modified_sites = [ind + 1 for ind, m in enumerate(site_iter) if "[80]" in m.group()]
                alt_sites = np.concatenate(ascore.alt_sites)
                n_overlapping_sites = np.intersect1d(modified_sites, alt_sites).shape[0]
                self.assertEqual(n_overlapping_sites, 0)

        test_consistency("velos_matches_1_mods.pkl", "velos_spectra_1_mods.pkl")
        test_consistency("velos_matches_2_mods.pkl", "velos_spectra_2_mods.pkl")
        test_consistency("velos_matches_3_mods.pkl", "velos_spectra_3_mods.pkl")

    def test_pepscore_return(self):
        ascore = PyAscore(bin_size=100., n_top=10,
                          mod_group="STY", mod_mass=79.966331)
        def test_return(match_file, spec_file):
            with open(os.path.join("test", "match_spectra_pairs", match_file), "rb") as src:
                match_list = pickle.load(src)

            with open(os.path.join("test", "match_spectra_pairs", spec_file), "rb") as src:
                spectra_list = pickle.load(src)

            for match, spectra in zip(match_list, spectra_list):
                ascore = PyAscore(bin_size=100., n_top=10,
                                  mod_group="STY", mod_mass=79.966331)
                ascore.score(spectra["mz_values"],
                             spectra["intensity_values"],
                             match["peptide"],
                             len(match["mod_positions"]))

                pep_scores = ascore.pep_scores
                nslots = len(pep_scores[0]["signature"])
                nmods = len(match["mod_positions"])
                ncombinations = np.math.factorial(nslots)/np.math.factorial(nmods)/np.math.factorial(nslots - nmods)
                self.assertEqual(
                    len(pep_scores), 
                    ncombinations
                )

                self.assertTrue(
                    np.all(np.diff([p["weighted_score"] for p in pep_scores]) <= 0)
                )

        test_return("velos_matches_1_mods.pkl", "velos_spectra_1_mods.pkl")
        test_return("velos_matches_2_mods.pkl", "velos_spectra_2_mods.pkl")
        test_return("velos_matches_3_mods.pkl", "velos_spectra_3_mods.pkl")

    def test_ambiguity(self):
        ascore = PyAscore(bin_size=100., n_top=10,
                          mod_group="STY", mod_mass=79.966331)
        def test_ambiguity(match_file, spec_file):
            with open(os.path.join("test", "match_spectra_pairs", match_file), "rb") as src:
                match_list = pickle.load(src)

            with open(os.path.join("test", "match_spectra_pairs", spec_file), "rb") as src:
                spectra_list = pickle.load(src)

            for match, spectra in zip(match_list, spectra_list):
                ascore = PyAscore(bin_size=100., n_top=10,
                                  mod_group="STY", mod_mass=79.966331)
                ascore.score(spectra["mz_values"],
                             spectra["intensity_values"],
                             match["peptide"],
                             len(match["mod_positions"]))

                pep_scores = ascore.pep_scores
                self.assertEqual(
                    ascore.calculate_ambiguity(pep_scores[0], pep_scores[0]), 0.
                )

                if len(pep_scores) > 1:
                    manual_ambiguity = ascore.calculate_ambiguity(pep_scores[0], pep_scores[1])
                    self.assertTrue(np.any(np.isclose(manual_ambiguity, ascore.ascores)))

        test_ambiguity("velos_matches_1_mods.pkl", "velos_spectra_1_mods.pkl")
        test_ambiguity("velos_matches_2_mods.pkl", "velos_spectra_2_mods.pkl")
        test_ambiguity("velos_matches_3_mods.pkl", "velos_spectra_3_mods.pkl")

