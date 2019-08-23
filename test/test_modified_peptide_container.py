import unittest
import numpy as np
from pyascore import PyModifiedPeptide

class TestPyModifiedPeptide(unittest.TestCase):
    def test_signature_incr(self):
        pep = PyModifiedPeptide("STY", 79.966331)

        # Only one sig available
        pep.consume_peptide("ASK", 1)
        self.assertFalse( pep.incr_signature() )

        pep.consume_peptide("PASSEFK", 2)
        self.assertFalse( pep.incr_signature() )

        # Small Sequence
        pep.consume_peptide("ASTK", 1)
        self.assertTrue( np.all( pep.get_signature() == np.array([1, 0]) ) )
        self.assertTrue( pep.incr_signature() ) # Incrementing should return 1 if a signature is available
        self.assertTrue( np.all( pep.get_signature() == np.array([0, 1]) ) )
        self.assertFalse( pep.incr_signature() ) # Incrementing should return 0 if no signature is available

        # Longer Sequence
        pep.consume_peptide("PASSSSSEFK", 2)
        self.assertTrue( np.all( pep.get_signature() == [1, 1, 0, 0, 0] ) )
        # increment several steps
        [pep.incr_signature() for i in range(4)]
        self.assertTrue( np.all( pep.get_signature() == [0, 1, 1, 0, 0] ) )
        
        # Make sure that incrementing stops as some point
        ind = 0
        while (pep.incr_signature() and ind < 1000):
            ind += 1
        self.assertTrue(ind < 1000)
        self.assertTrue( np.all( pep.get_signature() == [0, 0, 0, 0, 0] ) ) # All zeros is the end right now

    def test_fragment_incr(self):

        pep = PyModifiedPeptide("STY", 79.966331)
        
        # All fragments in first signature
        pep.consume_peptide("ASTK", 1)
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

        self.assertFalse(pep.incr_fragment())

        # Test correct fragment mz with signature increment
        pep.incr_signature()
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

        self.assertFalse(pep.incr_fragment())
        self.assertFalse(pep.incr_signature())
