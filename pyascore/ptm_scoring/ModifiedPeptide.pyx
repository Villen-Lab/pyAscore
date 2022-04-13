# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from libcpp.vector cimport vector
from libcpp.string cimport string
from ModifiedPeptide cimport ModifiedPeptide

cdef class PyModifiedPeptide:
    """
    The PyModifiedPeptide object provides functionality for modified residues of peptides.

    Objects can take in a sequence, a set of fixed position modifications, and a variable amount
    of unlocalized modifications which can fall on any residue of a specified type. The design 
    allows peaks from spectra to be matched to theoretical peaks from any possible localization. 
    Individual realizations of modified peptides are encoded via a "signature", which is merely a 
    binary vector with an entry for each modifiable residue, and a 1 signifying that the residue 
    is modified. The peaks of all possible modifications states can be traversed by creating a 
    PyFragmentGraph object, and if two signatures are provided, one can retrieve only the site 
    determining peaks.

    Parameters
    ----------
    mod_group : str
        A string which lists the possible modified residues for the unlocalized modification. For example, 
        with phosphorylation, you may want "STY".
    mod_mass : float
        The mass of the unlocalized modification in Daltons. For example, phosphorylation is 79.966331.
    mz_error : float
        The error in daltons to match theoretical peaks to consumed spectral peaks. The option to use PPM
        will likely be included in the future. (Defaults to 0.5)
    """
    cdef ModifiedPeptide * modified_peptide_ptr

    def __cinit__(self, str mod_group, float mod_mass, float mz_error = .5, str fragment_types = "by"):
        self.modified_peptide_ptr = new ModifiedPeptide(mod_group.encode("utf8"), 
                                                        mod_mass, 
                                                        mz_error, 
                                                        fragment_types.encode("utf8"))

    def __dealloc__(self):
        del self.modified_peptide_ptr

    def add_neutral_loss(self, str group, float mass):
        self.modified_peptide_ptr[0].addNeutralLoss(group.encode("utf8"), mass)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def consume_peptide(self, str peptide, size_t n_of_mod,
                        size_t max_fragment_charge=1,
                        np.ndarray[unsigned int, ndim=1, mode="c"] aux_mod_pos = None, 
                        np.ndarray[float, ndim=1, mode="c"] aux_mod_mass = None):
        """Consumes a single peptide sequence and creates it's internal representation

        Parameters
        ----------
        peptide : str
            The peptide string without any modifications or n-terminal markings
        n_of_mod : int > 0
            Number of unlocalized modifications on the sequence
        max_fragment_charge : int > 0
            Fragments will be considered from charge 1 to max_fragment_charge
        aux_mod_pos : ndarray of uint32
            Positions of fixed modifications. Most modification positions should start at 1 with 0 being
            reserved for n-terminal modifications, as seems to be the field prefered encoding.
        aux_mod_mass : ndarray of float32
            Masses of individual fixed postion modifications.
        """
        if aux_mod_pos is not None and aux_mod_mass is not None:
            self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"),
                                                        n_of_mod,
                                                        max_fragment_charge,
                                                        &aux_mod_pos[0], 
                                                        &aux_mod_mass[0], 
                                                        aux_mod_pos.size)
        else:
            self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod, max_fragment_charge)

    def get_peptide(self, np.ndarray[unsigned int, ndim=1, mode="c"] signature = np.ndarray(0, dtype=np.uint32)):
        """Prints the modified sequence (residues plus mod mass) of the consumed peptide

        Parameters
        ----------
        signature : ndarray of 0,1 values
            Encodes the modification state that each modifiable amino acid should have. Defaults to no modifications.

        Returns
        -------
        str
            A peptide sequence with bracketed modification masses, e.g. PEPT[80]IDEK.
        """
        cdef vector[size_t] signature_vector
        cdef size_t ind
        for ind in range(<size_t> signature.size):
            signature_vector.push_back(<size_t> signature[ind])

        return self.modified_peptide_ptr[0].getPeptide(signature_vector).decode("utf8")

    def get_fragment_graph(self, str fragment_type, size_t charge_state, str mode = "all"):
        """Builds a PyFragmentGraph object which references the current PyModifiedPeptide object

        Parameters
        ----------
        fragment_type : char
            The type of fragment graph to create, e.g. 'b'.
        charge_state : integer > 0
            The charge state of all fragments.

        Returns
        -------
        PyFragmentGraph
            The fragment graph of specified type and charge state.
        """
        return PyFragmentGraph(self, fragment_type.encode("utf8")[0], charge_state, mode)

    def get_site_determining_ions(self, np.ndarray[unsigned int, ndim=1, mode="c"] sig_1,
                                        np.ndarray[unsigned int, ndim=1, mode="c"] sig_2,
                                        str fragment_type, size_t max_charge):
        """Determine the non-overlapping theoretical fragments of two peptides.

        Parameters
        ----------
        sig1 : ndarray of 0,1 values
            Encodes the modification state that each modifiable amino acid should have in the first peptide.
        sig2 : ndarray of 0,1 values
            Encodes the modification state that each modifiable amino acid should have in the second peptide.
        fragment_type : char
            The type of fragment graph to create, e.g. 'b'.
        max_charge : integer > 0
            Site determining ions are produced for all charge states from 1 to max_charge inclusive.

        Returns
        -------
        tuple of ndarray
            A tuple of length 2 with entries that contain all theoretical fragments for a modified peptide not found
            in the other modified peptide.
        """

        cdef vector[size_t] sig_vec_1
        cdef vector[size_t] sig_vec_2
        cdef size_t ind
        for ind in range(<size_t> min(sig_1.size, sig_2.size)):
            sig_vec_1.push_back(<size_t> sig_1[ind])
            sig_vec_2.push_back(<size_t> sig_2[ind])

        cdef vector[vector[float]] ions = self.modified_peptide_ptr[0].getSiteDeterminingIons(
            sig_vec_1, sig_vec_2, fragment_type.encode("utf8")[0], max_charge
        )

        ion_arrays = (np.zeros(ions[0].size(), dtype=np.float32),
                      np.zeros(ions[1].size(), dtype=np.float32))
        for ind in range(ions[0].size()):
            ion_arrays[0][ind] = ions[0][ind]
        for ind in range(ions[1].size()):
            ion_arrays[1][ind] = ions[1][ind]

        return ion_arrays

cdef class PyFragmentGraph:
    """
    The PyFragmentGraph object allows traversal of the modification tree of a PyModifiedPeptide object.

    Every possible modified residue creates a branch point determined by it being modified or not, and 
    the PyFragmentGraph object allows efficient depth first traversal of the tree. By specifying whether
    the graph should be made over the b or y ions (more ion types to come), this object will spit out the
    appropriate theoretical MZ for each fragment. These MZ can be mathed against the internal cache of the 
    PyModifiedPeptide object to determine if any consumed peaks match the theoretical peak. For efficiency's
    sake, y and b ion graphs iterate through signatures differently, and may not necessarily be reverse 
    iterators of each other.

    Object creation through the PyModifiedPeptide.get_fragment_graph method is suggested by not required.

    Parameters
    ----------
    peptide : PyModifiedPeptide
        A PyModifiedPeptide instance which has consumed at least one peptide
    fragment_type : char
        The type of fragment graph to create, e.g. 'b'.
    charge_state : int > 0
        The charge state of all fragments

    Attributes
    ----------
    fragment_type : char
        The current type of fragment returned by the graph
    charge_state : int
        The current charge state for fragments returned from the graph
    """
    cdef str mode
    cdef ModifiedPeptide.FragmentGraph * fragment_graph_ptr

    def __cinit__(self, PyModifiedPeptide peptide, char fragment_type, size_t charge_state, str mode = "all"):
        assert mode in ("all", "reduced")
        self.mode = mode
        self.fragment_graph_ptr = new ModifiedPeptide.FragmentGraph(
            peptide.modified_peptide_ptr, fragment_type, charge_state
        )

    def __dealloc__(self):
        del self.fragment_graph_ptr

    @property
    def fragment_type(self):
        cdef bytes frag_type = self.fragment_graph_ptr[0].getFragmentType()
        return frag_type.decode('utf-8')

    @property
    def charge_state(self):
        return self.fragment_graph_ptr[0].getChargeState()

    def reset_iterator(self):
        """Resets iterator to the first position of the first signature."""
        return self.fragment_graph_ptr[0].resetIterator()

    def incr_signature(self):
        """Get next signature at position of last modification switch."""
        return self.fragment_graph_ptr[0].incrSignature()

    def is_signature_end(self):
        """Check if iterator is at last signature.

        Returns
        -------
        bool
            Is this the last signature?
        """
        return self.fragment_graph_ptr[0].isSignatureEnd()

    def reset_fragment(self):
        """Resets iterator to the first position of the current signature."""
        return self.fragment_graph_ptr[0].resetFragment()

    def incr_fragment(self):
        """Increment to next fragment for current signature."""
        return self.fragment_graph_ptr[0].incrFragment()

    def is_fragment_end(self):
        """Check if iterator has reached the last fragemnt, i.e. the end of the peptide."""
        return self.fragment_graph_ptr[0].isFragmentEnd()

    def set_signature(self, np.ndarray[unsigned int, ndim=1, mode="c"] new_signature):
        """Change signature to user specified value and reset to the first fragment.

        Parameters
        ----------
        new_signature : ndarray of uint32
            Encodes the modification state that each modifiable amino acid should have.
        """
        cdef vector[size_t] signature_vector
        cdef size_t ind
        for ind in range(<size_t> new_signature.size):
            signature_vector.push_back(<size_t> new_signature[ind])

        self.fragment_graph_ptr[0].setSignature(signature_vector)
        
    def get_signature(self):
        """Return current signature.

        Returns
        -------
        ndarray of uint64
            Array with one position per modifiable amino acid and a 1 if modified and 0 if not.
        """
        cdef vector[size_t] signature_vector = self.fragment_graph_ptr[0].getSignature()
        signature_array = np.zeros(signature_vector.size(), dtype=np.uint64)

        cdef size_t i = 0
        for i in range(signature_vector.size()):
            signature_array[i] = signature_vector[i]
        return signature_array

    def get_fragment_mz(self):
        """Return the size of the current fragment in m/z.

        Returns
        -------
        float
        """
        return self.fragment_graph_ptr[0].getFragmentMZ()

    def get_fragment_size(self):
        """Return the size of the current fragment in number of amino acids.

        Returns
        -------
        int
        """
        return self.fragment_graph_ptr[0].getFragmentSize()

    def get_fragment_seq(self):
        """Return sequence of current fragment without modifications.

        Returns
        -------
        str
        """
        return self.fragment_graph_ptr[0].getFragmentSeq().decode('utf8')

    def iter_permutations(self):
        """Iterate through remaining signatures and return fragment graph ready for iteration.

        If mode == 'all', reset fragments to the first position before returning graph.

        Yields
        ------
        PyFragmentGraph
            reference to current graph
        """
        while not self.is_signature_end():
            yield self
  
            self.incr_signature()
            if self.mode == "all":
              self.reset_fragment()

    def iter_fragments(self):
        """Iterate through remaining fragments of current signature.
        
        Yields
        ------
        (float, string)
            pair of fragment mz and fragment label
        """
        while not self.is_fragment_end():
            label = self.fragment_type + str(self.get_fragment_size())
            result = (self.get_fragment_mz(), label)
            self.incr_fragment()
            yield result
