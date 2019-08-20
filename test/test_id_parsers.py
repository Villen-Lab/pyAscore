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
