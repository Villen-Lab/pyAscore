# distutils: language = c++
import cython
import numpy as np
cimport numpy as np

from Ascore cimport Ascore, ScoreContainer
from ModifiedPeptide cimport ModifiedPeptide
from Spectra cimport BinnedSpectra

from libcpp.vector cimport vector

cdef class PyAscore:
    """
    The PyAscore object scores the localization of post translational modifications (PTMs).

    Objects are designed to take in a spectra, the associated peptide sequence, a set of fixed 
    position modifications, and a variable amount of unlocalized modifications and determine
    how much evidence exists for placing PTMs on individual amino acids. The algorithm is a 
    modified version of Beausoleil et al. [PMID: 16964243] which can efficiently handle any
    size peptide and arbitrary PTM masses. Each scored PSM will generate the most likely PTM
    positions and scores, as well as alternative sites for each PTM which have equal evidence
    but evidence that is less than or equal to the maximum. These alternative sites are not 
    required to be adjacent (i.e. not separated by another modifiable residue).

    Note:
        Attributes are only meaningful after consumption of the first peptide.

    Parameters
    ----------
    bin_size : float
        Size in MZ of each bin
    n_top : int
        Number of top peaks to retain in each bin (must be >= 0)
    mod_group : str
        A string which lists the possible modified residues for the unlocalized modification. For example, 
        with phosphorylation, you may want "STY".
    mod_mass : float
        The mass of the unlocalized modification in Daltons. For example, phosphorylation is 79.966331.
    mz_error : float
        The error in daltons to match theoretical peaks to consumed spectral peaks. The option to use PPM
        will likely be included in the future. (Defaults to 0.5)
    fragment_types : str
        The theoretical fragment types to score.

    Attributes
    ----------
    best_sequence : str
        Peptide sequence with modifications included in brackets for the best scoring localization.
    best_score : float
        The best PepScore among all possible localization permutations.
    pep_scores : list of dict
        Python dict representations of internal PepScore objects. Each object contains the sequence
        of the underlying peptide and all information necessary to calculate the ambiguity scores.
        The list is sorted by decreasing weighted_score which is also known as the PepScore.
    ascores : ndarray of float32
        Ascores for each individual non-static site in the peptide.
    alt_sites : list of ndarry of uint32
        Alternative positions for each individual non-static site in the peptide.
    """
    cdef Ascore * ascore_ptr
    cdef ModifiedPeptide * modified_peptide_ptr
    cdef BinnedSpectra * binned_spectra_ptr

    def __cinit__(self, float bin_size, size_t n_top,
                        str mod_group, float mod_mass,
                        float mz_error=.5,
                        str fragment_types="by"):
        self.binned_spectra_ptr = new BinnedSpectra(bin_size, n_top)
        self.modified_peptide_ptr = new ModifiedPeptide(mod_group.encode("utf8"), 
                                                        mod_mass, 
                                                        mz_error,
                                                        fragment_types.encode("utf8"))
        self.ascore_ptr = new Ascore()

    def __dealloc__(self):
        """Explicitly delete internal Cython pointers"""
        del self.binned_spectra_ptr
        del self.modified_peptide_ptr
        del self.ascore_ptr

    def add_neutral_loss(self, str group, float mass):
        """Add a neutral loss ion to any fragment containing specified amino acids

        For every fragment ion containing one of the amino acids defined in `group`,
        the algorithm will search for a secondary neutral loss peak to also score.
        Neutral losses can be specified for any amino acid using its one letter code,
        i.e STY, and can be placed on a modified amino acid (variable or otherwise)
        by using a lowercase letter, i.e. sty.

        Parameters
        ----------
        group : str
            Amino acids that should allow the neutral loss specified with their single
            letter code, i.e. STY for Ser, Thr, and Tyr.
        mass : float
            Absolute mass of the neutral loss. Negative values will cause the algorithm
            to look for higher mass peaks.
        """
        self.modified_peptide_ptr[0].addNeutralLoss(group.encode("utf8"), mass)

    @cython.boundscheck(False)
    @cython.wraparound(False)
    def score(self, np.ndarray[double, ndim=1, mode="c"] mz_arr not None, 
                    np.ndarray[double, ndim=1, mode="c"] int_arr not None,
                    str peptide, size_t n_of_mod,
                    size_t max_fragment_charge=1,
                    np.ndarray[unsigned int, ndim=1, mode="c"] aux_mod_pos = None,
                    np.ndarray[float, ndim=1, mode="c"] aux_mod_mass = None):
        """Consume spectra and associated peptide information and score PTM localization

        Parameters
        ----------
        mz_arr : ndarray of float64
            Array of MZ values for each peak in a spectra.
        int_arr : ndarray of float64
            Array of intensity values for each peak in a spectra.
        peptide : str
            The peptide string without any modifications or n-terminal markings.
        n_of_mod : int > 0
            Number of unlocalized modifications on the sequence.
        max_fragment_charge : int > 0
            Maximum fragment charge to be used for score calculations.
        aux_mod_pos : ndarray of uint32
            Positions of fixed modifications. Most modification positions should start at 1 with 0 being
            reserved for n-terminal modifications, as seems to be the field prefered encoding.
        aux_mod_mass : ndarray of float32
            Masses of individual fixed postion modifications.
        """ 
        # Consume spectra and bin
        self.binned_spectra_ptr[0].consumeSpectra(&mz_arr[0], &int_arr[0], mz_arr.size)

        # Build modified peptide with or without constant mods
        if aux_mod_pos is not None and aux_mod_mass is not None:
            self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod,
                                                        max_fragment_charge,
                                                        &aux_mod_pos[0], &aux_mod_mass[0],
                                                        aux_mod_pos.size)
        else:
            self.modified_peptide_ptr[0].consumePeptide(peptide.encode("utf8"), n_of_mod, max_fragment_charge)
        
        # Allow modified peptide to consume peaks from binned spectra
        while (self.binned_spectra_ptr[0].getBin() < self.binned_spectra_ptr[0].getNBins()):

            self.binned_spectra_ptr[0].resetRank()
            while (self.binned_spectra_ptr[0].getRank() < self.binned_spectra_ptr[0].getNPeaks()):
                self.modified_peptide_ptr[0].consumePeak(self.binned_spectra_ptr[0].getMZ(),
                                                         self.binned_spectra_ptr[0].getRank())
                self.binned_spectra_ptr[0].nextRank()

            self.binned_spectra_ptr[0].nextBin()

        self.ascore_ptr[0].score(self.binned_spectra_ptr[0], self.modified_peptide_ptr[0])

    def _score_cont_to_pyobj(self, ScoreContainer cont):
        """Internal function to convert C++ ScoreContainer to native python dict
        
        Returns
        -------
        dict
            Python dictionary with attributes: signature, counts, scores, weighted_score, total_fragments
        """
        cdef size_t i
        pyobj = {}

        pyobj["signature"] = np.zeros(cont.signature.size(),
                                      dtype=np.int32)
        for i in range(cont.signature.size()):
            pyobj["signature"][i] = cont.signature[i]
 
        pyobj["counts"] = np.zeros(cont.counts.size(),
                                   dtype=np.int32)
        pyobj["scores"] = np.zeros(cont.scores.size(),
                                   dtype=np.float32)
        for i in range(cont.counts.size()):
            pyobj["counts"][i] = cont.counts[i]
            pyobj["scores"][i] = cont.scores[i]
 
        pyobj["weighted_score"] = cont.weighted_score
        pyobj["total_fragments"] = cont.total_fragments

        return pyobj

    def _pyobj_to_score_cont(self, dict pyobj):
        """Internal function to convert native python dict to C++ ScoreContainer
        
        Returns
        -------
        ScoreContainer
            C++ object which can be transfered to internal algorithms
        """
        cdef size_t i, nitems
        cdef ScoreContainer cont

        nitems = pyobj["signature"].shape[0]
        for i in range(nitems):
            cont.signature.push_back(pyobj["signature"][i])

        nitems = pyobj["counts"].shape[0]
        for i in range(nitems):
            cont.counts.push_back(pyobj["counts"][i])
            cont.scores.push_back(pyobj["scores"][i])

        cont.weighted_score = pyobj["weighted_score"]
        cont.total_fragments = pyobj["total_fragments"]

        return cont

    def calculate_ambiguity(self, dict ref_score, dict other_score):
        """Calculate ambiguity between 2 competing localizations

        Inputs to this function should come directly from the pep_scores attribute.
        For this score to be possitive, the score dict with the highest weighted_score
        should come first. When the weighted_score for both is equal, this will return 0.

        Parameters
        ----------
        ref_score : dict
            PepScore object for one localization.
        other_score : dict
            PepScore object for competing localization.

        Returns
        -------
        float
            Relative evidence for localization in first argument vs the second argument
        """
        cdef ScoreContainer ref_score_cont = self._pyobj_to_score_cont(ref_score)
        cdef ScoreContainer other_score_cont = self._pyobj_to_score_cont(other_score)

        return self.ascore_ptr[0].calculateAmbiguity(ref_score_cont, other_score_cont)

    @property
    def best_sequence(self):
        return self.ascore_ptr[0].getBestSequence().decode("utf8")

    @property
    def best_score(self):
        return self.ascore_ptr[0].getBestScore()

    @property
    def pep_scores(self):
        cdef vector[ScoreContainer] raw_score_conts = self.ascore_ptr[0].getAllPepScores();
        cdef vector[string] sequences = self.ascore_ptr[0].getAllSequences();

        proc_score_conts = []
        cdef size_t i
        for i in range(raw_score_conts.size()):
            pyobj_cont = self._score_cont_to_pyobj(raw_score_conts[i])
            pyobj_cont["sequence"] = sequences[i].decode("utf8")
            proc_score_conts.append(pyobj_cont)

        return proc_score_conts;

    @property
    def ascores(self):
        cdef vector[float] score_vector = self.ascore_ptr[0].getAscores()
        cdef np.ndarray[float, ndim=1, mode="c"] score_array = np.zeros(
            score_vector.size(), dtype=np.float32
        )

        cdef size_t i = 0
        for i in range(score_vector.size()):
            score_array[i] = score_vector[i]
        return score_array

    @property
    def alt_sites(self):
        cdef size_t nmods
        cdef size_t mod_ind
        cdef size_t alt_ind
        
        cdef vector[size_t] alt_vector
        cdef np.ndarray[np.uint32_t, ndim=1, mode="c"] alt_array

        alt_site_list = []

        n_mods = self.modified_peptide_ptr[0].getNumberOfMods()
        for mod_ind in range(n_mods):
            alt_vector = self.ascore_ptr[0].getAlternativeSites(mod_ind)

            alt_array = np.zeros( alt_vector.size(), dtype=np.uint32 )
            
            for alt_ind in range(alt_vector.size()):
                alt_array[alt_ind] = alt_vector[alt_ind]

            alt_site_list.append(alt_array)

        return alt_site_list
