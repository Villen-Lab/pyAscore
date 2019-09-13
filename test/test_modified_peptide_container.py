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

    def test_fragment_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)
        
        def test(graph, charge, neutral_masses):
            for m in neutral_masses:
                m = (m + charge) / max(1, charge)
                self.assertTrue(
                    np.isclose(graph.get_fragment_mz(), m, rtol=1e-6, atol=0)
                )
                graph.incr_fragment()

        pep.consume_peptide("ASTK", 1)
        max_charge = 3
        for c in range(max_charge + 1):
            b_graph = pep.get_fragment_graph("b", c)
            # First signature all the way through
            test(b_graph, c, np.array([71.03711 , 238.035471, 339.083151, 467.178111]))
            self.assertTrue(b_graph.is_fragment_end())
            # Second signature, picking up from common node
            b_graph.incr_signature()
            test(b_graph, c, np.array([158.06914 , 339.083151, 467.178111]))
            # Second signature, from beginning
            b_graph.reset_iterator()
            b_graph.incr_signature()
            test(b_graph, c, np.array([71.03711, 158.06914 , 339.083151, 467.178111]))

            y_graph = pep.get_fragment_graph("y", c)
            # First signature all the way through
            test(y_graph, c, np.array([146.11024 , 327.124251, 414.156281, 485.193391]))
            self.assertTrue(y_graph.is_fragment_end())
            # Second signature, picking up from common node
            y_graph.incr_signature()          
            test(y_graph, c, np.array([247.15792 , 414.156281, 485.193391]))
            # Second signature, from beginning
            y_graph.reset_iterator()
            y_graph.incr_signature()
            test(y_graph, c, np.array([146.11024, 247.15792 , 414.156281, 485.193391]))
