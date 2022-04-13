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

        pep.consume_peptide("ASK", 1, 1, np.array([0]).astype(np.uint32), np.array([20.]).astype(np.float32))
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

        def test(graph, true_masses):
            for m in true_masses:
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        # Test b fragments
        b_graph = pep.get_fragment_graph("b", 1)
        b_graph.set_signature(np.array([0, 1, 1, 0, 0], dtype=np.uint32))
        true_b_masses = [98.06058, 169.09769, 256.12972, 423.12808, 590.12644,
                         677.15847, 764.19050, 893.23309, 1040.30150]
        test(b_graph, true_b_masses)

        # Test c fragments
        c_graph = pep.get_fragment_graph("c", 1)
        c_graph.set_signature(np.array([0, 1, 1, 0, 0], dtype=np.uint32))
        true_c_masses = [115.08713, 186.12424, 273.15627, 440.15463, 607.15299,
                         694.18502, 781.21705, 910.25964, 1057.32805]
        test(c_graph, true_c_masses)

        # Test y fragments
        y_graph = pep.get_fragment_graph("y", 1)
        y_graph.set_signature(np.array([0, 1, 1, 0, 0], dtype=np.uint32))
        true_y_masses = [147.11334, 294.18176, 423.22435, 510.25638, 597.28841,
                         764.28677, 931.28513, 1018.3171, 1089.3542]
        test(y_graph, true_y_masses)

        # Test z fragments
        z_graph = pep.get_fragment_graph("z", 1)
        z_graph.set_signature(np.array([0, 1, 1, 0, 0], dtype=np.uint32))
        true_z_masses = [130.08680, 277.15521, 406.19780, 493.22983, 580.26186,
                         747.26022, 914.25858, 1001.29061, 1072.32772]
        test(z_graph, true_z_masses)

    def test_iterator(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        # Test mode == "all"
        pep.consume_peptide("ASMTK", 1)
        b_graph = pep.get_fragment_graph("b", 1)
        true_sig_list = [[1, 0], [0, 1]]
        true_frag_list = [np.array([71.03711 , 238.03547, 369.07596, 470.12364]),
                          np.array([71.03711, 158.06914 , 289.10963, 470.12364])]
        for graph, true_sig, true_frags in zip(b_graph.iter_permutations(), true_sig_list, true_frag_list):
            self.assertTrue(np.all(graph.get_signature() == true_sig))
            for (pred_frag, label), mass in zip(graph.iter_fragments(), true_frags):
                self.assertTrue(
                    np.isclose(pred_frag, mass + 1.007825, rtol=1e-6, atol=0)
                )

        # Test mode == "reduced"
        pep.consume_peptide("ASMTK", 1)
        b_graph = pep.get_fragment_graph("b", 1, mode="reduced")
        true_sig_list = [[1, 0], [0, 1]]
        true_frag_list = [np.array([71.03711 , 238.03547, 369.07596, 470.12364]),
                          np.array([158.06914 , 289.10963, 470.12364])]
        for graph, true_sig, true_frags in zip(b_graph.iter_permutations(), true_sig_list, true_frag_list):
            self.assertTrue(np.all(graph.get_signature() == true_sig))
            for (pred_frag, label), mass in zip(graph.iter_fragments(), true_frags):
                self.assertTrue(
                    np.isclose(pred_frag, mass + 1.007825, rtol=1e-6, atol=0)
                )

    def test_fragment_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)
        
        def test(graph, charge, neutral_masses):
            for m in neutral_masses:
                m = (m + charge * 1.007825) / max(1, charge)
                #print(graph.get_fragment_mz(), m)
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        # Test unmodified peptide
        pep.consume_peptide("ASMTK", 1)
        max_charge = 3
        for c in range(max_charge + 1):
            # Test b fragments
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([71.03711 , 238.03547, 369.07596, 470.12364]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([158.06914 , 289.10963, 470.12364]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([71.03711, 158.06914 , 289.10963, 470.12364]))
            b_graph.reset_fragment()
            test(b_graph, c, np.array([71.03711, 158.06914 , 289.10963, 470.12364]))

            # Test c fragments
            c_graph = pep.get_fragment_graph("c", c)
            # First signature all the way through
            test(c_graph, c, np.array([88.06365, 255.06201, 386.10251, 487.15019]))
            self.assertTrue(c_graph.is_fragment_end())
            # Second signature, picking up from common node
            c_graph.incr_signature()
            test(c_graph, c, np.array([175.09568, 306.13617, 487.15019]))
            # Second signature, from beginning
            c_graph.reset_iterator()
            c_graph.incr_signature()
            test(c_graph, c, np.array([88.06365, 175.09568, 306.13617, 487.15019]))
            c_graph.reset_fragment()
            test(c_graph, c, np.array([88.06365, 175.09568, 306.13617, 487.15019]))

            # Test y fragments
            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.10552, 327.11953, 458.160025, 545.19205]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()          
            test(y_graph, c, np.array([247.15320, 378.19369, 545.192056]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.105525, 247.15320, 378.19369, 545.192056]))
            y_graph.reset_fragment()
            test(y_graph, c, np.array([146.105525, 247.15320, 378.19369, 545.192056]))

            # Test z fragments
            z_graph = pep.get_fragment_graph("z", c)
            # First signature all the way through
            test(z_graph, c, np.array([129.07897, 310.09298, 441.13347, 528.16550]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            z_graph.incr_signature()
            test(z_graph, c, np.array([230.12665, 361.16714, 528.16550]))
            # Second signature, from beginning
            z_graph.reset_iterator()
            z_graph.incr_signature()
            test(z_graph, c, np.array([129.07897, 230.12665, 361.16714, 528.16550]))
            z_graph.reset_fragment()
            test(z_graph, c, np.array([129.07897, 230.12665, 361.16714, 528.16550]))

            # Test Z fragments
            z_graph = pep.get_fragment_graph("Z", c)
            # First signature all the way through
            test(z_graph, c, np.array([130.086795, 311.100805, 442.141295, 529.173325]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            z_graph.incr_signature()
            test(z_graph, c, np.array([231.134475, 362.174965, 529.173325]))
            # Second signature, from beginning
            z_graph.reset_iterator()
            z_graph.incr_signature()
            test(z_graph, c, np.array([130.086795, 231.134475, 362.174965, 529.173325]))
            z_graph.reset_fragment()
            test(z_graph, c, np.array([130.086795, 231.134475, 362.174965, 529.173325]))

    def test_fragment_incr_terminal(self):
        pep = PyModifiedPeptide("nKc", 42.010565)

        def test(graph, charge, neutral_masses):
            for m in neutral_masses:
                m = (m + charge * 1.007825) / max(1, charge)
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        # Test unmodified peptide
        pep.consume_peptide("ASKTR", 1)
        max_charge = 3
        for c in range(max_charge + 1):
            # Test b fragments
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([113.047675, 200.079705, 328.174664, 429.222344]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([71.03711, 158.06914, 328.174664, 429.222344]))
            # Third signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([286.16409, 387.21178]))

            # Test y fragments
            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([216.12223, 317.16992, 445.26487, 532.29691]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()
            test(y_graph, c, np.array([174.11167, 275.15935, 445.26487, 532.2969]))
            # Third signature, picking up from common node
            y_graph.incr_signature()
            test(y_graph, c, np.array([403.254314, 490.286345]))

    def test_fragment_incr_mod(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        def test(graph, charge, neutral_masses):
            for m in neutral_masses:
                m = (m + charge * 1.007825) / max(1, charge)
                #print(graph.get_fragment_mz(), m)
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        # Test modified peptide
        pep.consume_peptide("ASMTK", 1, 1, np.array([0, 3]).astype(np.uint32), np.array([42.010565, 15.994915]).astype(np.float32))
        max_charge = 3
        for c in range(max_charge + 1):
            # Test b fragments
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([113.04767, 280.04603, 427.08144, 528.12912]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([200.07970, 347.11510, 528.12912]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([113.04767, 200.07970, 347.11510, 528.12912]))

            # Test c fragments
            c_graph = pep.get_fragment_graph("c", c)
            # First signature all the way through
            test(c_graph, c, np.array([130.07422, 297.07258, 444.10798, 545.15567]))
            self.assertTrue(c_graph.is_fragment_end())
            # Second signature, picking up from common node
            c_graph.incr_signature()
            test(c_graph, c, np.array([217.10625, 364.14165, 545.15567]))
            # Second signature, from beginning
            c_graph.reset_iterator()
            c_graph.incr_signature()
            test(c_graph, c, np.array([130.07422, 217.10625, 364.14165, 545.15567]))
         
            # Test y fragments
            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.10552, 327.11953, 474.15494, 561.18697]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()
            test(y_graph, c, np.array([247.15320, 394.18861, 561.18697]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.10552, 247.15320, 394.18861, 561.18697]))

            # Test z fragments
            z_graph = pep.get_fragment_graph("z", c)
            # First signature all the way through
            test(z_graph, c, np.array([129.07897, 310.09298, 457.12839, 544.16042]))
            self.assertTrue(z_graph.is_fragment_end())
            # Second signature, picking up from common node
            z_graph.incr_signature()
            test(z_graph, c, np.array([230.12665, 377.16206, 544.16042]))
            # Second signature, from beginning
            z_graph.reset_iterator()
            z_graph.incr_signature()
            test(z_graph, c, np.array([129.07897, 230.12665, 377.16206, 544.16042]))

            # Test Z fragments
            z_graph = pep.get_fragment_graph("Z", c)
            # First signature all the way through
            test(z_graph, c, np.array([130.086795, 311.100805, 458.136215, 545.168245]))
            self.assertTrue(z_graph.is_fragment_end())
            # Second signature, picking up from common node
            z_graph.incr_signature()
            test(z_graph, c, np.array([231.134475, 378.169885, 545.168245]))
            # Second signature, from beginning
            z_graph.reset_iterator()
            z_graph.incr_signature()
            test(z_graph, c, np.array([130.086795, 231.134475, 378.169885, 545.168245]))

    def test_fragment_incr_nl(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        def test(graph, charge, neutral_masses):
            for m in neutral_masses:
                m = (m + charge * 1.007825) / max(1, charge)
                #print(graph.get_fragment_mz(), m)
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        pep.add_neutral_loss("ST", 18.01528)
        pep.consume_peptide("ASMTK", 1)
        max_charge = 3
        for c in range(max_charge + 1):
            # Test b fragments
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([71.03711, 238.035471, 369.075961, 
                                       470.123641, 452.108361])) 
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([158.06914, 140.05386, 289.10963, 
                                       271.09435, 470.123641, 452.108361]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([71.03711, 158.06914, 140.05386, 289.10963, 
                                       271.09435, 470.123641, 452.108361]))

            # Test y fragments
            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.10552, 327.11953, 458.16002, 
                                       545.19205, 527.17677]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()
            test(y_graph, c, np.array([247.15320, 229.13792, 378.19369, 360.17841, 545.19205, 527.17677]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.10552, 247.15320, 229.13792, 378.19369, 
                                       360.17841, 545.19205, 527.17677]))

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
        true_mz = (np.array([234.14537, 365.18586]),
                   np.array([314.11171, 445.15220]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        # Site determining peaks with a single straddling constant modification
        pep.consume_peptide("ASMSK", 1, 1,
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
        true_mz = (np.array([234.14537, 381.18078]),
                   np.array([314.11171, 461.14711]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

    def test_site_determining_double(self):
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
        true_mz = (np.array([982.35929]),                                          
                   np.array([1062.32562]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        # Site determining peaks with one interfering mod
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
        true_mz = (np.array([510.25638, 597.28841, 728.32890, 895.32726, 982.35929]),
                   np.array([590.22271, 677.25474, 808.29523, 975.29359, 1062.32562])) 
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        # Site determining peaks with two interfering mods
        pep.consume_peptide("PASSSMSSEFK", 2, 1,
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
        true_mz = (np.array([510.25638, 597.28841, 744.32381, 911.32217, 998.35420]), 
                   np.array([590.22271, 677.25474, 824.29014, 991.28850, 1078.32053]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

    def test_site_determining_high_charge(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        pep.consume_peptide("ASMHSK", 1, 2)
        calc_mz = pep.get_site_determining_ions(np.array([1, 0], dtype=np.uint32),
                                                np.array([0, 1], dtype=np.uint32),
                                                "b", 2)
        true_mz = (np.array([120.02556, 185.545805, 239.0427475, 254.07526, 370.083786, 507.142696]),
                   np.array([80.042395, 145.56264, 159.076965, 214.092095, 290.117455, 427.176365]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

        calc_mz = pep.get_site_determining_ions(np.array([1, 0], dtype=np.uint32),
                                                np.array([0, 1], dtype=np.uint32),
                                                "y", 2)
        true_mz = (np.array([117.576602, 186.106057, 234.14537, 251.626302, 371.20428, 502.24477]),
                   np.array([157.559767, 226.089222, 291.609468, 314.111710, 451.17062, 582.211111]))
        for c, m in zip(calc_mz, true_mz):
            self.assertTrue(np.all(np.isclose(c, m, rtol=1e-05, atol=0)))

    def test_peptide_print(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        pep.consume_peptide("ASMTK", 1, 1, np.array([0, 3]).astype(np.uint32), np.array([42.010565, 15.994915]).astype(np.float32))
        self.assertEqual(pep.get_peptide(), "n[42]AS[80]M[16]TK")
        self.assertEqual(pep.get_peptide(np.array([0, 1], dtype=np.uint32)), "n[42]ASM[16]T[80]K")

        pep.consume_peptide("PASSSSSEFK", 2)
        self.assertEqual(pep.get_peptide(), "PAS[80]S[80]SSSEFK")
        self.assertEqual(pep.get_peptide(np.array([0, 1, 0, 1, 0], dtype=np.uint32)), "PASS[80]SS[80]SEFK")
