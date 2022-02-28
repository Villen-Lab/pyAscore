import unittest
import os
import pickle
from itertools import product
from pyascore import id_parsers
import numpy as np
from pyteomics import mass

STD_AA_MASS = mass.std_aa_mass

class TestMassCorrector(unittest.TestCase):
    corrector = id_parsers.MassCorrector()

    def test_n_term(self):
        res = "X"
        mass = 42.010565
        for i in range(6):
            
            c_res, c_pos, c_mass = self.corrector.correct(res, 0, round(mass, i))
            self.assertEqual(
                (c_res[0], c_pos[0], c_mass[0]),
                ('n', 0, mass)
            )

        res = "M"
        n_mod_mass = 42.010565
        mass = STD_AA_MASS[res] + n_mod_mass
        for i in range(6):

            c_res, c_pos, c_mass = self.corrector.correct(res, 1, round(mass, i))
            self.assertEqual(
                (c_res[0], c_pos[0], c_mass[0]),
                ('n', 0, n_mod_mass)
            )

    def test_n_term_combined(self):
        res = "M"
        n_mod_mass = 42.010565
        oxi_mass = 15.9949
        mass = STD_AA_MASS[res] + n_mod_mass + oxi_mass
        for i in range(6):

            c_res, c_pos, c_mass = self.corrector.correct(res, 1, round(mass, i))
            self.assertEqual(
                (c_res[0], c_pos[0], c_mass[0]),
                ('n', 0, n_mod_mass)
            )
            
            self.assertEqual(
                (c_res[1], c_pos[1], c_mass[1]),
                ('M', 1, oxi_mass)
            )

    def test_res(self):
        res = "S"
        phospho_mass = 79.966331
        mass = STD_AA_MASS[res] + phospho_mass
        for i in range(6):

            c_res, c_pos, c_mass = self.corrector.correct(res, 5, round(mass, i))
            self.assertEqual(
                (c_res[0], c_pos[0], c_mass[0]),
                (res, 5, phospho_mass)
            )

    def test_not_found(self):
        res = "M"
        phospho_mass = 79.966331
        mass = STD_AA_MASS[res] + phospho_mass
        for i in range(6):

            try:
                c_res, c_pos, c_mass = self.corrector.correct(res, 5, round(mass, i))

            except ValueError:
                continue

    def test_multiple(self):
        n_mod_mass = 42.010565
        oxi_mass = 15.9949
        phospho_mass = 79.966331

        peptide = "MRAMSLVSNEGDSEQNEIR"
        uncorrected_positions = np.array([1, 5])
        uncorrected_masses = np.array([STD_AA_MASS["M"] + n_mod_mass + oxi_mass,
                                       STD_AA_MASS["S"] + phospho_mass])

        true_positions = np.array([0, 1, 5])
        true_masses = np.array([n_mod_mass,
                                oxi_mass,
                                phospho_mass])

        corrected_positions, corrected_masses = self.corrector.correct_multiple(peptide,
                                                                                uncorrected_positions,
                                                                                uncorrected_masses)

        self.assertTrue(np.all(corrected_positions == true_positions),
                        "Positions are {}, not {}".format(corrected_positions, true_positions))
        self.assertTrue(np.all(corrected_masses == true_masses),
                        "Masses are {}, not {}".format(corrected_positions, true_positions))


def example_generator(file_name):
    with open(file_name, "rb") as source:
        examples = pickle.load(source)
        for e in examples:
            yield e


class TestIDExtractors(unittest.TestCase):

    program_list = ["comet", "percolator"]
    instrument_list = ["qexactive", "velos"]
    

    global_answers = {("comet", "qexactive") : [
                          dict(scan=2, charge_states=2, peptides="MRAMSLVSNEGDSEQNEIR", mod_positions=np.array([ 1,  5,  8, 13])),
                          dict(scan=3, charge_states=2, peptides="KEESEESDDDMGFGLFD", mod_positions=np.array([ 4, 7, 11 ])),
                          dict(scan=4, charge_states=2, peptides="KEESEESDDDMGFGLFD", mod_positions=np.array([ 4, 7 ]))
                      ],
                      ("comet", "velos") : [
                          dict(scan=2, charge_states=3, peptides="QADIQSTVLQINMPRGDLPVGNYQKMAKLADAR", mod_positions=np.array([ 13, 23 ])),
                          dict(scan=3, charge_states=4, peptides="ALSTCASHFTAVSVFYGTVIFIYLQPSSSHSMDTDK", mod_positions=np.array([ 5, 10, 28, 32 ])),
                          dict(scan=4, charge_states=2, peptides="LLVKKIVSLVR", mod_positions=np.array([]))
                      ],
                      ("percolator", "qexactive") : [
                          dict(scan=26840, charge_states=3, peptides="ATVPVAAATAAEGEGSPPAVAAVAGPPAAAEVGGGVGGSSR", mod_positions=np.array([ 16 ])),
                          dict(scan=27795, charge_states=2, peptides="GEADLFDSGDIFSTGTGSQSVER", mod_positions=np.array([ 16 ])),
                          dict(scan=22462, charge_states=3, peptides="LAEAPSPAPTPSPTPVEDLGPQTSTSPGR", mod_positions=np.array([]))
                      ],
                      ("percolator", "velos") : [
                          dict(scan=28126, charge_states=3, peptides="KGDVVHCWYTGTLQDGTVFDTNIQTSAK", mod_positions=np.array([ 7 ])),
                          dict(scan=33362, charge_states=3, peptides="HQILEQAVEDYAETVHQLSK", mod_positions=np.array([])),
                          dict(scan=28509, charge_states=3, peptides="RMATEVAADALGEEWKGYVVR", mod_positions=np.array([]))
                      ],
                     }

    def test_pepxml_extractor(self):
        extractor = id_parsers.PepXMLExtractor()

        for prog, instr in product(self.program_list, self.instrument_list):
            file_name = "_".join([prog, instr, "pepxml", "examples"]) + ".pkl"

            for ind, examp in enumerate(example_generator(
                                  os.path.join("test", "pyteomics_examples", "pepxml", file_name)
                              )):
                extracted_data = extractor.extract(examp)
                answers = self.global_answers[(prog, instr)]
                
                self.assertEqual(extracted_data["scans"][0], answers[ind]["scan"])
                self.assertEqual(extracted_data["charge_states"][0], answers[ind]["charge_states"])
                self.assertEqual(extracted_data["peptides"][0], answers[ind]["peptides"])
                self.assertTrue(np.all(extracted_data["peptides"][0] == answers[ind]["peptides"])) # Comparing arrays


class TestIdnetificationParser(unittest.TestCase):

    tide_answers = [{"scan": 14760, "charge_state": 3, "score": 4.48925829, "peptide": "KMSDDEDDDEEEYGKEEHEK",
                     "mod_positions": np.array([3]), "mod_masses": np.array([79.966331])},
                    {"scan": 14760, "charge_state": 3, "score": 1.60863817, "peptide": "KMSDDEDDDEEEYGKEEHEK",
                     "mod_positions": np.array([13]), "mod_masses": np.array([79.966331])},
                    {"scan": 18330, "charge_state": 3, "score": 2.81006455, "peptide": "EDLPAENGETKTEESPASDEAGEK",
                     "mod_positions": np.array([18]), "mod_masses": np.array([79.966331])},
                    {"scan": 18330, "charge_state": 3, "score": 2.54836297, "peptide": "EDLPAENGETKTEESPASDEAGEK",
                     "mod_positions": np.array([15]), "mod_masses": np.array([79.966331])},
                    {"scan": 20462, "charge_state": 3, "score": 3.30343485, "peptide": "RRASWASENGETDAEGTQMTPAK",
                     "mod_positions": np.array([4]), "mod_masses": np.array([79.966331])},
                    {"scan": 20462, "charge_state": 3, "score": 2.74408746, "peptide": "RRASWASENGETDAEGTQMTPAK",
                     "mod_positions": np.array([7]), "mod_masses": np.array([79.966331])},
                    {"scan": 21996, "charge_state": 3, "score": 3.88421941, "peptide": "AEEPPSQLDQDTQVQDMDEGSDDEEEGQK",
                     "mod_positions": np.array([17, 21]), "mod_masses": np.array([15.9949  , 79.966331])},
                    {"scan": 21996, "charge_state": 3, "score": 2.35982895, "peptide": "AEEPPSQLDQDTQVQDMDEGSDDEEEGQK",
                     "mod_positions": np.array([12, 17]), "mod_masses": np.array([79.966331, 15.9949  ])},
                    {"scan": 26219, "charge_state": 3, "score": 3.62030768, "peptide": "GKEELAEAEIIKDSPDSPEPPNK",
                     "mod_positions": np.array([17]), "mod_masses": np.array([79.966331])},
                    {"scan": 26219, "charge_state": 3, "score": 3.5756743, "peptide": "GKEELAEAEIIKDSPDSPEPPNK",
                     "mod_positions": np.array([14]), "mod_masses": np.array([79.966331])},
                    {"scan": 26962, "charge_state": 3, "score": 3.6686945, "peptide": "KEDSDEEEDDDSEEDEEDDEDEDEDEDEIEPAAMK",
                     "mod_positions": np.array([ 4, 12]), "mod_masses": np.array([79.966331, 79.966331])},
                    {"scan": 26962, "charge_state": 3, "score": 0.7472344, "peptide": "NASNVKHHDSSALGVYSYIPLVENPYFSSWPPSGTSSK",
                     "mod_positions": np.array([11, 29]), "mod_masses": np.array([79.966331, 79.966331])},
                    {"scan": 27845, "charge_state": 3, "score": 2.49656582, "peptide": "DLGSTEDGDGTDDFLTDKEDEK",
                     "mod_positions": np.array([16]), "mod_masses": np.array([79.966331])},
                    {"scan": 27845, "charge_state": 3, "score": 2.15698767, "peptide": "DLGSTEDGDGTDDFLTDKEDEK",
                     "mod_positions": np.array([11]), "mod_masses": np.array([79.966331])},
                    {"scan": 31328, "charge_state": 3, "score": 5.46427345, "peptide": "EGHSLEMENENLVENGADSDEDDNSFLK",
                     "mod_positions": np.array([ 7, 19]), "mod_masses": np.array([15.9949  , 79.966331])},
                    {"scan": 31328, "charge_state": 3, "score": 4.58486271, "peptide": "EGHSLEMENENLVENGADSDEDDNSFLK",
                     "mod_positions": np.array([ 7, 25]), "mod_masses": np.array([15.9949  , 79.966331])},
                    {"scan": 32257, "charge_state": 3, "score": 3.3405838, "peptide": "KPATPAEDDEDDDIDLFGSDNEEEDK",
                     "mod_positions": np.array([ 4, 19]), "mod_masses": np.array([79.966331, 79.966331])},
                    {"scan": 32257, "charge_state": 3, "score": 1.01620567, "peptide": "AGDMGNCVSGQQQEGGVSEEMKGPVQEDK",
                     "mod_positions": np.array([7]), "mod_masses": np.array([57.021464])},
                    {"scan": 35669, "charge_state": 3, "score": 3.87117219, "peptide": "VEEESTGDPFGFDSDDESLPVSSK",
                     "mod_positions": np.array([14]), "mod_masses": np.array([79.966331])},
                    {"scan": 35669, "charge_state": 3, "score": 3.24795341, "peptide": "VEEESTGDPFGFDSDDESLPVSSK",
                     "mod_positions": np.array([18]), "mod_masses": np.array([79.966331])}]

    def test_pepxml_reader(self):
        target_file="test/example_inputs/psms/test_psms.pep.xml"
        parser = id_parsers.IdentificationParser(target_file, "pepXML", score_string="xcorr_score")
        psm_list = parser.to_list()

        self.assertEqual(len(psm_list), 20)
        for psm, answer in zip(psm_list, self.tide_answers):
            self.assertEqual(psm["scan"], answer["scan"])
            self.assertEqual(psm["charge_state"], answer["charge_state"])
            self.assertEqual(psm["score"], answer["score"])
            self.assertEqual(psm["peptide"], answer["peptide"])
            self.assertTrue(np.all(psm["mod_positions"] == answer["mod_positions"]))
            self.assertTrue(np.all(psm["mod_masses"] == answer["mod_masses"]))
   
    def test_mzidml_reader(self):
        target_file="test/example_inputs/psms/test_psms.mzid"
        parser = id_parsers.IdentificationParser(target_file, "mzIdentML", score_string="SEQUEST:xcorr")
        psm_list = parser.to_list()

        self.assertEqual(len(psm_list), 20)
        # mzIdentMLs from tide collapse localizations,
        # so the following will only look at every other hit.
        for psm, answer in zip(psm_list[::2], self.tide_answers[::2]):
            self.assertEqual(psm["scan"], answer["scan"])
            self.assertEqual(psm["charge_state"], answer["charge_state"])
            self.assertEqual(psm["score"], answer["score"])
            self.assertEqual(psm["peptide"], answer["peptide"])
            self.assertTrue(np.all(psm["mod_positions"] == answer["mod_positions"]))
            self.assertTrue(np.all(psm["mod_masses"] == answer["mod_masses"]))
