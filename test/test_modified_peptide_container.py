import unittest
import numpy as np
from itertools import product
from pyascore import PyModifiedPeptide

class TestPyModifiedPeptide(unittest.TestCase):
    def test_one_sig_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        def test(pep):
            graph_b = pep.get_fragment_graph("b", 1)
            graph_b.incr_signature()
            self.assertTrue(graph_b.is_signature_end())

            graph_y = pep.get_fragment_graph("y", 1)
            graph_y.incr_signature()
            self.assertTrue(graph_y.is_signature_end())

        pep.consume_peptide("ASK", 1)
        test(pep)

        pep.consume_peptide("PASSEFK", 2)
        test(pep)

        pep.consume_peptide("ASK", 1, np.array([0]).astype(np.uint32), np.array([20.]).astype(np.float32))
        test(pep)

    def test_signature_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        # Small Sequence
        pep.consume_peptide("ASTK", 1)

        graph_b_progress = [np.array([1, 0]), np.array([0, 1])]
        graph_y_progress = [np.array([0, 1]), np.array([1, 0])]
        graph_b = pep.get_fragment_graph("b", 1)
        graph_y = pep.get_fragment_graph("y", 1)
        for ind in range(2):
            self.assertTrue( np.all( graph_b.get_signature() == graph_b_progress[ind] ) )
            self.assertTrue( np.all( graph_y.get_signature() == graph_y_progress[ind] ) )
            graph_b.incr_signature(), graph_y.incr_signature()
        self.assertTrue(graph_b.is_signature_end() and graph_y.is_signature_end())

        # Longer Sequence
        pep.consume_peptide("PASSSSSEFK", 2)

        graph_b_progress = [np.array([1, 1, 0, 0, 0]), np.array([0, 1, 1, 0, 0])]
        graph_y_progress = [np.array([0, 0, 0, 1, 1]), np.array([0, 0, 1, 1, 0])]
        graph_b = pep.get_fragment_graph("b", 1)
        graph_y = pep.get_fragment_graph("y", 1)
        for ind in range(2):
            self.assertTrue( np.all( graph_b.get_signature() == graph_b_progress[ind] ) )
            self.assertTrue( np.all( graph_y.get_signature() == graph_y_progress[ind] ) )
            [ (graph_b.incr_signature(), graph_y.incr_signature()) for x in range(4)]
        self.assertFalse(graph_b.is_signature_end() or graph_y.is_signature_end())

    def test_signature_stop(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        pep.consume_peptide("PASSSSSEFK", 2)
        graph_b = pep.get_fragment_graph("b", 1)
        graph_y = pep.get_fragment_graph("y", 1)
        while (not graph_b.is_signature_end() or not graph_y.is_signature_end()):
            graph_b.incr_signature(), graph_y.incr_signature()
        self.assertTrue( np.all( graph_b.get_signature() == [0, 0, 0, 0, 0] ) )
        self.assertTrue( np.all( graph_y.get_signature() == [0, 0, 0, 0, 0] ) )

    def test_set_signature(self):
        seq = "PASSSSSEFK"
        pep = PyModifiedPeptide("STY", 79.966331)
        pep.consume_peptide("PASSSSSEFK", 2)

        b_graph = pep.get_fragment_graph("b", 1)
        b_graph.set_signature(np.array([0, 1, 1, 0, 0], dtype=np.uint32))
        true_b_masses = [97.05276, 168.08987, 255.1219, 422.120261, 589.118622, 676.150652, 
                         763.182682, 892.225272, 1039.293682, 1167.388642]
        for m in true_b_masses:
            m = (m + 1.00727647)
            self.assertTrue(
                np.isclose(b_graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
            )
            b_graph.incr_fragment()


        y_graph = pep.get_fragment_graph("y", 1)
        y_graph.set_signature(np.array([0, 1, 1, 0, 0], dtype=np.uint32))
        true_y_masses = [146.11024, 293.17865, 422.22124, 509.25327, 596.2853, 
                         763.283661, 930.282022, 1017.314052, 1088.351162, 1185.403922]
        for m in true_y_masses:
            m = (m + 1.00727647)
            self.assertTrue(
                np.isclose(y_graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
            )
            y_graph.incr_fragment()

    def test_fragment_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)
        
        def test(graph, charge, neutral_masses):
            for m in neutral_masses:
                m = (m + charge * 1.00727647) / max(1, charge)
                #print(graph.get_fragment_mz(), m)
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        pep.consume_peptide("ASMTK", 1)
        max_charge = 3
        for c in range(max_charge + 1):
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([71.03711 , 238.035471, 369.075961, 470.123641, 598.218601]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([158.06914 , 289.10963, 470.123641, 598.218601]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([71.03711, 158.06914 , 289.10963, 470.123641, 598.218601]))

            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.11024, 327.124251, 458.164741, 545.196771, 616.233881]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()          
            test(y_graph, c, np.array([247.15792, 378.19841, 545.196771, 616.233881]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.11024, 247.15792, 378.19841, 545.196771, 616.233881]))

        pep.consume_peptide("ASMTK", 1, np.array([0, 3]).astype(np.uint32), np.array([42.010565, 15.994915]).astype(np.float32))
        max_charge = 3
        for c in range(max_charge + 1):
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([113.047675, 280.046036, 427.081441, 528.129121, 656.224081]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([200.079705, 347.11510999999996, 528.1291209999999, 656.224081]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([113.047675, 200.079705, 347.11510999999996, 528.1291209999999, 656.224081]))

            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.11024, 327.124251, 474.159656, 561.191686, 674.239361]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()
            test(y_graph, c, np.array([247.15792, 394.193325, 561.191686, 674.239361]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.11024, 247.15792, 394.193325, 561.191686, 674.239361]))

        pep.consume_peptide("ASMECTK", 1, np.array([3, 5]).astype(np.uint32), np.array([15.994915, 57.021464]).astype(np.float32))
        max_charge = 3
        for c in range(max_charge + 1):
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([71.03711, 238.035471, 385.070876, 514.113466, 674.14412, 775.1918, 903.28676]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([158.06914, 305.104545, 434.147135, 594.177789, 775.1918, 903.28676]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([71.03711, 158.06914, 305.104545, 434.147135, 594.177789, 775.1918, 903.28676]))

            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.11024, 327.124251, 487.154905, 616.197495, 763.2329, 850.26493, 921.30204]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()
            test(y_graph, c, np.array([247.15792, 407.1885734, 536.231164, 683.266569, 850.26493, 921.30204]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.11024, 247.15792, 407.1885734, 536.231164, 683.266569, 850.26493, 921.30204]))

    def test_site_determining_single(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        # Site determining peaks with no constant modifications
        pep.consume_peptide("ASMSK", 1)
        calc_mz = pep.get_site_determining_ions(np.array([1, 0], dtype=np.uint32),
                                                np.array([0, 1], dtype=np.uint32),
                                                "b", 1)
        true_mz = (np.array([239.0427475, 370.08323747]), 
                   np.array([159.07641647, 290.11690647]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        calc_mz = pep.get_site_determining_ions(np.array([1, 0], dtype=np.uint32),
                                                np.array([0, 1], dtype=np.uint32),
                                                "y", 1)
        true_mz = (np.array([234.14954647, 365.19003647]),
                   np.array([314.11587747, 445.15636747]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        # Site determining peaks with a single straddling constant modification
        pep.consume_peptide("ASMSK", 1,
                            np.array([3], dtype=np.uint32),
                            np.array([15.9949146202], dtype=np.float32))
        calc_mz = pep.get_site_determining_ions(np.array([1, 0], dtype=np.uint32),
                                                np.array([0, 1], dtype=np.uint32),
                                                "b", 1)
        true_mz = (np.array([239.0427475, 386.078152]), 
                   np.array([159.07641647, 306.111821]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        calc_mz = pep.get_site_determining_ions(np.array([1, 0], dtype=np.uint32),
                                                np.array([0, 1], dtype=np.uint32),
                                                "y", 1)
        true_mz = (np.array([234.14954647, 381.184951]),
                   np.array([314.11587747, 461.1512821]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

    def test_sit_determining_double(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        # Site determining peaks with no interfering mods
        pep.consume_peptide("PASSSMSSEFK", 2)
        calc_mz = pep.get_site_determining_ions(np.array([1, 0, 0, 1, 0], dtype=np.uint32),
                                                np.array([0, 1, 0, 1, 0], dtype=np.uint32),
                                                "b", 1)
        true_mz = (np.array([336.09550747]), 
                   np.array([256.12917647]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        calc_mz = pep.get_site_determining_ions(np.array([1, 0, 0, 1, 0], dtype=np.uint32),
                                                np.array([0, 1, 0, 1, 0], dtype=np.uint32),
                                                "y", 1)
        true_mz = (np.array([982.36345747]),                                          
                   np.array([1062.32978847]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        # Site determining peaks with one interfering mods
        pep.consume_peptide("PASSSMSSEFK", 2)
        calc_mz = pep.get_site_determining_ions(np.array([1, 0, 1, 0, 0], dtype=np.uint32),
                                                np.array([0, 0, 1, 0, 1], dtype=np.uint32),
                                                "b", 1)
        true_mz = (np.array([336.09550747, 423.12753747, 590.12589847, 721.16638847, 808.19841847]),                                          
                   np.array([256.12917647, 343.16120647, 510.15956747, 641.20005747, 728.23208747]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        calc_mz = pep.get_site_determining_ions(np.array([1, 0, 1, 0, 0], dtype=np.uint32),
                                                np.array([0, 0, 1, 0, 1], dtype=np.uint32),
                                                "y", 1)
        true_mz = (np.array([510.26054647, 597.29257647, 728.33306647, 895.33142747, 982.36345747]), 
                   np.array([590.22687747, 677.25890747, 808.29939747, 975.29775847, 1062.32978847]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        # Site determining peaks with two interfering mods
        pep.consume_peptide("PASSSMSSEFK", 2,
                            np.array([6], dtype=np.uint32),
                            np.array([15.9949146202], dtype=np.float32))
        calc_mz = pep.get_site_determining_ions(np.array([1, 0, 1, 0, 0], dtype=np.uint32),
                                                np.array([0, 0, 1, 0, 1], dtype=np.uint32),
                                                "b", 1)
        true_mz = (np.array([336.09550747, 423.12753747, 590.12589847, 737.1613030902, 824.1933330902]),
                   np.array([256.12917647, 343.16120647, 510.15956747, 657.1949720902, 744.2270020902]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        calc_mz = pep.get_site_determining_ions(np.array([1, 0, 1, 0, 0], dtype=np.uint32),
                                                np.array([0, 0, 1, 0, 1], dtype=np.uint32),
                                                "y", 1)
        true_mz = (np.array([510.26054647, 597.29257647, 744.3279810902, 911.32634209, 998.35837209]),  
                   np.array([590.22687747, 677.25890747, 824.2943121, 991.29267309, 1078.324703]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

    def test_peptide_print(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        pep.consume_peptide("ASMTK", 1, np.array([0, 3]).astype(np.uint32), np.array([42.010565, 15.994915]).astype(np.float32))
        self.assertEqual(pep.get_peptide(), "n[42]AS[80]M[16]TK")
        self.assertEqual(pep.get_peptide(np.array([0, 1], dtype=np.uint32)), "n[42]ASM[16]T[80]K")

        pep.consume_peptide("PASSSSSEFK", 2)
        self.assertEqual(pep.get_peptide(), "PAS[80]S[80]SSSEFK")
        self.assertEqual(pep.get_peptide(np.array([0, 1, 0, 1, 0], dtype=np.uint32)), "PASS[80]SS[80]SEFK")
