import unittest
import numpy as np
from pyascore import PyModifiedPeptide

class TestPyModifiedPeptide(unittest.TestCase):
    def test_signature_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        # Only one sig available
        pep.consume_peptide("ASK", 1)
        self.assertFalse( pep.incr_signature() )

        graph_b = pep.get_fragment_graph("b", 1)
        graph_b.incr_signature()
        self.assertTrue(graph_b.is_signature_end())

        graph_y = pep.get_fragment_graph("y", 1)
        graph_y.incr_signature()
        self.assertTrue(graph_y.is_signature_end())

        pep.consume_peptide("PASSEFK", 2)
        self.assertFalse( pep.incr_signature() )

        graph_b = pep.get_fragment_graph("b", 1)
        graph_b.incr_signature()
        self.assertTrue(graph_b.is_signature_end())

        graph_y = pep.get_fragment_graph("y", 1)
        graph_y.incr_signature()
        self.assertTrue(graph_y.is_signature_end())

        # Small Sequence
        pep.consume_peptide("ASTK", 1)
        self.assertTrue( np.all( pep.get_signature() == np.array([1, 0]) ) )
        self.assertTrue( pep.incr_signature() ) # Incrementing should return 1 if a signature is available
        self.assertTrue( np.all( pep.get_signature() == np.array([0, 1]) ) )
        self.assertFalse( pep.incr_signature() ) # Incrementing should return 0 if no signature is available

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
        self.assertTrue( np.all( pep.get_signature() == [1, 1, 0, 0, 0] ) )
        # increment several steps
        [pep.incr_signature() for i in range(4)]
        self.assertTrue( np.all( pep.get_signature() == [0, 1, 1, 0, 0] ) )

        graph_b_progress = [np.array([1, 1, 0, 0, 0]), np.array([0, 1, 1, 0, 0])]
        graph_y_progress = [np.array([0, 0, 0, 1, 1]), np.array([0, 0, 1, 1, 0])]
        graph_b = pep.get_fragment_graph("b", 1)
        graph_y = pep.get_fragment_graph("y", 1)
        for ind in range(2):
            self.assertTrue( np.all( graph_b.get_signature() == graph_b_progress[ind] ) )
            self.assertTrue( np.all( graph_y.get_signature() == graph_y_progress[ind] ) )
            [ (graph_b.incr_signature(), graph_y.incr_signature()) for x in range(4)]
        self.assertFalse(graph_b.is_signature_end() or graph_y.is_signature_end())
        
        # Make sure that incrementing stops as some point
        ind = 0
        while (pep.incr_signature() and ind < 1000):
            ind += 1
        self.assertTrue(ind < 1000)
        self.assertTrue( np.all( pep.get_signature() == [0, 0, 0, 0, 0] ) ) # All zeros is the end right now

        while (not graph_b.is_signature_end() or not graph_y.is_signature_end()):
            graph_b.incr_signature(), graph_y.incr_signature()
        self.assertTrue( np.all( graph_b.get_signature() == [0, 0, 0, 0, 0] ) )
        self.assertTrue( np.all( graph_y.get_signature() == [0, 0, 0, 0, 0] ) )

    def test_fragment_incr(self):

        pep = PyModifiedPeptide("STY", 79.966331)
        
        # All fragments in first signature
        max_charge = 3
        pep.consume_peptide("ASTK", 1)
        b_fragment_graphs = [pep.get_fragment_graph("b", charge) for charge in range(max_charge + 1)]
        y_fragment_graphs = [pep.get_fragment_graph("y", charge) for charge in range(max_charge + 1)]
        b_neutral_mass = np.array([71.03711 , 238.035471, 339.083151, 467.178111])
        y_neutral_mass = np.array([146.11024 , 327.124251, 414.156281, 485.193391])
        charge_1_mass = [72.03711, 239.035471, 340.083151, 468.17811099999994]
        charge_2_mass = [36.518555, 120.0177355, 170.5415755, 234.5890555]
        for i in range(4):

            self.assertTrue( 
                np.isclose(pep.get_fragment_mz(1), charge_1_mass[i], rtol=1e-6, atol=0)
            )
            self.assertTrue(
                np.isclose(pep.get_fragment_mz(2), charge_2_mass[i], rtol=1e-6, atol=0)
            )

            pep.incr_fragment()
 
            for charge in range(max_charge + 1):
                b_mass = (b_neutral_mass[i] + charge) / max(charge, 1)
                self.assertTrue( 
                    np.isclose( b_fragment_graphs[charge].get_fragment_mz(), b_mass, rtol=1e-7, atol=0 )
                )
                b_fragment_graphs[charge].incr_fragment()

                y_mass = (y_neutral_mass[i] + charge) / max(charge, 1)
                self.assertTrue(
                    np.isclose( y_fragment_graphs[charge].get_fragment_mz(), y_mass, rtol=1e-7, atol=0 )
                )
                y_fragment_graphs[charge].incr_fragment()

        self.assertFalse(pep.incr_fragment())
        # Test correct fragment mz with signature increment
        pep.incr_signature()
        [g.incr_signature() for g in b_fragment_graphs]
        [g.incr_signature() for g in y_fragment_graphs]
        b_neutral_mass = np.array([158.06914 , 339.083151, 467.178111])
        y_neutral_mass = np.array([247.15792 , 414.156281, 485.193391])
        charge_1_mass = [159.06914, 340.083151, 468.178111]
        charge_2_mass = [80.03457, 170.5415755, 234.5890555]
        for i in range(3):

            self.assertTrue(
                np.isclose(pep.get_fragment_mz(1), charge_1_mass[i], rtol=1e-6, atol=0) 
            )
            self.assertTrue(
                np.isclose(pep.get_fragment_mz(2), charge_2_mass[i], rtol=1e-6, atol=0)
            )

            pep.incr_fragment()

            for charge in range(max_charge + 1):
                b_mass = (b_neutral_mass[i] + charge) / max(charge, 1)
                self.assertTrue(
                    np.isclose( b_fragment_graphs[charge].get_fragment_mz(), b_mass, rtol=1e-7, atol=0 )
                )
                b_fragment_graphs[charge].incr_fragment()
                
                y_mass = (y_neutral_mass[i] + charge) / max(charge, 1)
                self.assertTrue(
                    np.isclose( y_fragment_graphs[charge].get_fragment_mz(), y_mass, rtol=1e-7, atol=0 )
                )
                y_fragment_graphs[charge].incr_fragment()

        self.assertFalse(pep.incr_fragment())
        self.assertFalse(pep.incr_signature())
