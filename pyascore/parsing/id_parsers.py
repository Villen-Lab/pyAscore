import re
import warnings
import numpy as np
from numpy import isclose
from pandas import read_csv
from pyteomics.mzid import MzIdentML
from pyteomics.pepxml import PepXML
from pyteomics import mass

# Constants
STD_AA_MASS = mass.std_aa_mass

# In case a set of custom mods is not supplied
COMMON_MODS = {"n": 42.010565,
               "M": 15.9949,
               "K": 8.014199,
               "S": 79.966331,
               "T": 79.966331,
               "Y": 79.966331,
               "C": 57.021464}

class MassCorrector:
    """
    A class to provide modification mass correction

    Often the mass of a modification is rounded when it comes from PSM files,
    or the mass of the modification may be combined with the mass of the amino
    acid and needs to be decouple. These situations make the mod difficult to
    work with, so this class corrects rounded and combined modification masses.

    Parameters
    ----------
    mod_mass_dict : dict
        A dictionary of known modifications to be used for correction
    aa_mass_dict : dict
        Correct masses for the amino acids in the peptide sequence
    mz_tol : float
        How far off can a rounded mass be to a known mod before it is counted as unknown
    n_mod_ind : int
        Is an n-terminal mod considered to be before the first AA (0) or after the first AA (1)
    """
    def __init__(self, mod_mass_dict=COMMON_MODS, aa_mass_dict=STD_AA_MASS, mz_tol=1.5, n_mod_ind=0):
        self.mod_mass_dict = mod_mass_dict
        self.aa_mass_dict = aa_mass_dict
        self.mz_tol = mz_tol
        self.n_mod_ind = n_mod_ind

    def correct(self, res, pos, mass):
        """Perform modification correction

        If the modification occurs on the n-terminus, and the n-term mod has been combined
        with another modification, this function will return tuples of length 2. Otherwise,
        a tuple of length 1 with be returned.

        Parameters
        ----------
        res : str
            The AA that the modification is on
        pos : int
            The position of the modification in the peptide
        mass : float
            The mass of the modification to be corrected

        Returns
        (str, ?str)
            A tuple of length 1 or 2 with the residues of the mod
        (int, ?int)
            A tuple of length 1 or 2 with the final positions of the mod
        (float, ?float)
            A tuple of length 1 or 2 with the final masses of the mod
        """
        # Store the modification masses if the can be found
        # otherwise, store inf, which is impossible to match
        std_mass = STD_AA_MASS.get(res, np.inf)
        mod_mass = self.mod_mass_dict.get(res, np.inf)
        n_mod_mass = self.mod_mass_dict.get('n', np.inf)

        # Tests various modification situations that seem
        # to occure depending on how various authors deem
        # it appropropriate to encode masses.
        if pos == 0 and isclose(mass, n_mod_mass,
                                 rtol=0., atol=self.mz_tol):
            return ('n',), (self.n_mod_ind,), (n_mod_mass,)

        elif pos == 1 and isclose(mass, std_mass + n_mod_mass,
                                  rtol=0., atol=self.mz_tol):
            return ('n',), (self.n_mod_ind,), (n_mod_mass,)

        elif pos == 1 and isclose(mass, std_mass + mod_mass + n_mod_mass,
                                  rtol=0., atol=self.mz_tol):
            return ('n', res), (self.n_mod_ind, pos), (n_mod_mass, mod_mass)

        elif isclose(mass, std_mass + mod_mass,
                     rtol=0., atol=self.mz_tol):
            return (res,), (pos,), (mod_mass,)

        else:
            pred_mod_mass = mass - STD_AA_MASS.get(res, 0.)
            warnings.warn("Unrecognized mod on {} at position {} with mass: {}"
                          " Using uncorrected mass.".format(res, pos, pred_mod_mass))
            return (res,), (pos,), (pred_mod_mass,)

    def correct_multiple(self, peptide, positions, masses):
        """Run the correction function on a list of mods for a single peptide

        Parameters
        ----------
        peptide : str
             The unmodified peptide sequence to reference
        positions : list
             A list of modification positions
        masses : list
             A list of modification masses

        Returns
        -------
        np.array
            An array of corrected positions
        np.array
            An array of corrected masses
        """
        corrected_positions = []
        corrected_masses = []
        for ind in range(len(positions)):
            pos = positions[ind]
            mass = masses[ind]
            res = 'n' if pos == 0 else peptide[pos - 1]
            _, pos, mass = self.correct(res, pos, mass)
            corrected_positions.extend(pos)
            corrected_masses.extend(mass)

        return np.array(corrected_positions), np.array(corrected_masses)

    def correct_numpy(self, peptide, positions, masses):
        """Correct masses taking advantage of numpy internals. DEPRECATED"""
        # Trivial return
        if positions.size == 0:
            return np.array([]), np.array([])

        corrected_positions = []
        corrected_masses = []

        n_mod_mass = self.mod_mass_dict.get('n', np.inf)
        n_mod_possibilities = [n_mod_mass,
                               STD_AA_MASS.get(peptide[0], np.inf) + n_mod_mass,
                               STD_AA_MASS.get(peptide[0], np.inf) + self.mod_mass_dict.get(peptide[0], np.inf) + n_mod_mass]
        if positions[0] == 0 and isclose(masses[0], n_mod_possibilities[0],
                                         rtol=0., atol=self.mz_tol):
            corrected_positions.append(0.)
            corrected_masses.append(n_mod_mass)
            positions = positions[1:]
            masses = masses[1:]

        elif positions[0] == 1 and isclose(
            masses[0], n_mod_possibilities[1],
            rtol=0., atol=self.mz_tol):
            corrected_positions.append(0.)
            corrected_masses.append(n_mod_mass)
            positions = positions[1:]
            masses = masses[1:]

        elif positions[0] == 1 and isclose(
            masses[0], n_mod_possibilities[2],
            rtol=0., atol=self.mz_tol):
            corrected_positions.append(0)
            corrected_masses.append(n_mod_mass)
            masses[0] -= n_mod_mass

        if positions.size > 0:
            std_masses = np.array([STD_AA_MASS.get(peptide[pos-1], np.inf) for pos in positions if pos != 0])
            std_mods = np.array([self.mod_mass_dict.get(peptide[pos-1], np.inf) for pos in positions if pos != 0])
            matches = np.isclose(masses, std_masses + std_mods, rtol=0., atol=self.mz_tol)
            if np.any(~matches):
                raise ValueError("Unrecognized mod at positions, {},"
                                 " with masses, {}".format(positions[np.where(~matches)],
                                                           masses[np.where(~matches)])) 
            corrected_positions.extend(positions)
            corrected_masses.extend(std_mods)    

        return np.array(corrected_positions), np.array(corrected_masses)


class PercolatorTXT:
    """
    Simple wrapper to provide standardized access to Percolator ouput csvs

    The only reason this is necessary is that I need access to this csv data to
    have the same interface as file access in pyteomics.

    Parameters
    ----------
    file_path : str
        Path to file to read
    """
    def __init__(self, file_path):
        self._internal_csv = read_csv(file_path, sep="\t")

    def map(self, func):
        """Map function over scans"""
        return self._internal_csv.groupby("scan").apply(func).to_list()


class MokapotTXT:
    """
    Simple wrapper to provide standardized access to Mokapot ouput csvs

    The only reason this is necessary is that I need access to this csv data to
    have the same interface as file access in pyteomics.

    Parameters
    ----------
    file_path : str
        Path to file to read
    """
    def __init__(self, file_path):
        self._internal_csv = read_csv(file_path, sep="\t")

    def map(self, func):
        """Map function over scans"""
        return self._internal_csv.groupby("ScanNr").apply(func).to_list()


class IDExtractor:
    """
    Parser class for peptide spectrum match objects

    This class defines the core public functionality of identification extraction.
    Users will pass single entry dictionaries via the extract method, which will
    then call parsing methods specific to each file type, before sanitizing final
    results relavant to modifications, and passing these back in a common form.

    Paramters:
    ----------
    score_string : str
        Specific string describing the score to extract from a record.
    static_mods : dict
        Dictionary of static mods to be used by some subclasses.
    """
    def __init__(self, score_string=None, static_mods={}):
        self._match_list = None
        self._match = None
        self.score_string = score_string
        self.entry = None
        self.static_mods = static_mods

    def _get_score(self):
        """Returns values requested via score_string if present
        
        Returns
        -------
        Score for individual PSM
        """
        try:
            return self._match[self.score_string]
        except KeyError:
            return None
    
    def _initialize_results(self):
        """Initializes the results dictionary with the correct sized list"""
        nmatches = len(self._match_list)
        fields = ["scans", "scores", "charge_states",
                  "peptides", "mod_positions", "mod_masses"]
        self.results = {f: [None] * nmatches for f in fields}

    def extract(self, entry):
        """Sanitizes a single entry by calling filetype specific extraction methods

        Parameters:
        -----------
        entry : dict
            A single dictionary like object of an identification file entry.
        """
        self.entry = entry
        self._match_list = self._get_matches()

        self._initialize_results()
        for ind, self._match in enumerate(self._match_list):
            self.results["scans"][ind] = self._get_scans()
            self.results["scores"][ind] = self._get_score()
            self.results["charge_states"][ind] = self._get_charge()
            self.results["peptides"][ind] = self._get_peptide()
            self.results["mod_positions"][ind], self.results["mod_masses"][ind] = self._get_mod_info()
        return self.results


class MzIdentMLExtractor(IDExtractor):
    def _get_matches(self):
        """Extract PSMs for a given scan
        
        Return
        ------
        list
            Matches for a given scan
        """
        try:
            return self.entry["SpectrumIdentificationItem"]
        except KeyError:
            return []

    def _get_scans(self):
        """Extract scan number

        Return
        ------
        int
            Scan number for PSM
        """
        try:
            # Specific search
            scan_search = re.search("(?<=scan=)([0-9]+)", self.entry["spectrumID"])

            # If specific search fails, do liberal first number search
            if scan_search is None:
                scan_search = re.search("[0-9]+", self.entry["spectrumID"])

            return int(scan_search.group())
        except KeyError:
            return -1

    def _get_charge(self):
        """Extract precursor charge

        Return
        ------
        int
           Charge of precursor
        """
        try:
            return int(self._match["chargeState"])
        except KeyError:
            return 0

    def _get_score(self):
        """Extract score based on score string

        Return
        ------
        float
            Score requested by user
        """
        try:
            return float(self._match[self.score_string])
        except KeyError:
            return None

    def _get_peptide(self):
        """Extract unmodified peptide sequence

        Returns
        -------
        str
            Peptide sequence without modifications
        """
        try:
            return self._match["PeptideSequence"]
        except KeyError:
            return ""

    def _get_mod_info(self):
        """Extract information about peptide modifications
        
        Return
        ------
        (np.array, np.array)
            Tuple with positions of modifications (0) and the masses of the modified AAs (1)
        """
        try:
            mod_list = self._match["Modification"]
            pos = np.zeros(len(mod_list), dtype=np.int32)
            mass = np.zeros(len(mod_list), dtype=np.float32)
            for ind, mod in enumerate(mod_list):
                pos[ind] = int(mod["location"])
                
                # Sometimes the amino acid is already available, sometimes not
                if "residues" in mod:
                    aa = mod["residues"][0]
                else:
                    aa = ("n" + self._get_peptide() + "c")[pos[ind]]
                
                mass[ind] = STD_AA_MASS.get(aa, 0.) + float(mod["monoisotopicMassDelta"])

            return pos, mass

        except KeyError:
            return np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.float32)


class PepXMLExtractor(IDExtractor):
    def _get_matches(self):
        """Extract PSMs for a given scan

        Return
        ------
        list
            Matches for a given scan
        """
        try:
            return self.entry["search_hit"]
        except KeyError:
            return []

    def _get_scans(self):
        """Extract scan number

        Return
        ------
        int
            Scan number for PSM
        """
        try:
            return int(self.entry["start_scan"])
        except KeyError:
            return -1

    def _get_charge(self):
        """Extract precursor charge

        Return
        ------
        int
           Charge of precursor
        """
        try:
            return int(self.entry["assumed_charge"])
        except KeyError:
            return 0

    def _get_score(self):
        """Extract score based on score string

        Return
        ------
        float
            Score requested by user
        """
        try:
            return float(self._match["search_score"][self.score_string])
        except KeyError:
            return None

    def _get_peptide(self):
        """Extract unmodified peptide sequence

        Returns
        -------
        str
            Peptide sequence without modifications
        """
        try:
            return self._match["peptide"]
        except KeyError:
            return ""

    def _get_mod_info(self):
        """Extract information about peptide modifications

        Return
        ------
        (np.array, np.array)
            Tuple with positions of modifications (0) and the masses of the modified AAs (1)
        """
        try:
            mod_list = self._match["modifications"]
            pos = np.zeros(len(mod_list), dtype=np.int32)
            mass = np.zeros(len(mod_list), dtype=np.float32)
            for ind, mod in enumerate(mod_list):
                pos[ind] = int(mod["position"])
                mass[ind] = float(mod["mass"])
            return pos, mass

        except KeyError:
            return np.zeros(0, dtype=np.int32), np.zeros(0, dtype=np.float32)


class PercolatorTXTExtractor(IDExtractor):
    def _get_matches(self):
        """Extract PSMs for a given scan

        This will transform pandas Series or Dataframe in list of dictionaries.

        Return
        ------
        list
            Matches for a given scan
        """
        return [id for id in self.entry.transpose().to_dict().values()]

    def _get_scans(self):
        """Extract scan number

        Return
        ------
        int
            Scan number for PSM
        """
        return self._match["scan"]

    def _get_charge(self):
        """Extract precursor charge

        Return
        ------
        int
           Charge of precursor
        """
        return self._match["charge"]

    def _get_score(self):
        """Extract percolator score

        Ignores score string inputed by user.

        Return
        ------
        float
            Percolator score
        """
        return self._match["percolator score"]

    def _get_peptide(self):
        """Extract unmodified peptide sequence

        Returns
        -------
        str
            Peptide sequence without modifications
        """
        return re.sub("n|(\[[^A-z]+\])", "", self._match["sequence"])

    def _get_mod_info(self):
        """Extract information about peptide modifications

        Return
        ------
        (np.array, np.array)
            Tuple with positions of modifications (0) and the masses of the modified AAs (1)
        """
        mass = []
        pos = []
        # Extract n terminus first
        n_term_mod = re.match("n?\[[^A-z]+\]", self._match["sequence"])
        if n_term_mod is not None:
            mass.append(
                re.search("(?<=\[)[^A-z]+(?=\])", n_term_mod.group()).group()
            )
            pos.append(0)
        elif "n" in self.static_mods:
            mass.append(self.static_mods["n"])
            pos.append(0)

        # Extract the rest of the sequence
        residues = [res.group() for res in re.finditer("[A-Z](\[[^A-z]+\])?", self._match["sequence"])]
        for ind, res in enumerate(residues, 1):
            extracted_mod = re.search("(?<=\[)[^A-z]+(?=\])", res)
            if extracted_mod is not None:
                mass.append(STD_AA_MASS[res[0]] + float(extracted_mod.group()))
                pos.append(ind)
            elif res[0] in self.static_mods:
                mass.append(STD_AA_MASS[res[0]] + self.static_mods[res[0]])
                pos.append(ind)
        return np.array(pos, dtype=np.int32), np.array(mass, dtype=np.float32)


class MokapotTXTExtractor(IDExtractor):
    def _get_matches(self):
        """Extract PSMs for a given scan

        This will transform pandas Series or Dataframe in list of dictionaries.

        Return
        ------
        list
            Matches for a given scan
        """
        return [id for id in self.entry.transpose().to_dict().values()]

    def _get_scans(self):
        """Extract scan number

        Return
        ------
        int
            Scan number for PSM
        """
        return self._match["ScanNr"]

    def _get_charge(self):
        """Extract precursor charge

        Return
        ------
        int
           Charge of precursor
        """
        return None

    def _get_score(self):
        """Extract mokapot score

        Ignores score string inputed by user.

        Return
        ------
        float
            Mokapot score
        """
        return self._match["mokapot score"]

    def _get_peptide(self):
        """Extract unmodified peptide sequence

        Returns
        -------
        str
            Peptide sequence without modifications
        """
        seq = re.sub("(^.\.)|(\..$)", "", self._match["Peptide"])
        return re.sub("n|(\[[^A-z]+\])", "", seq)

    def _get_mod_info(self):
        """Extract information about peptide modifications

        Return
        ------
        (np.array, np.array)
            Tuple with positions of modifications (0) and the masses of the modified AAs (1)
        """
        mass = []
        pos = []
        seq = re.sub("(^.\.)|(\..$)", "", self._match["Peptide"])
        # Extract n terminus first
        n_term_mod = re.match("n?\[[^A-z]+\]", seq)
        if n_term_mod is not None:
            mass.append(
                re.search("(?<=\[)[^A-z]+(?=\])", n_term_mod.group()).group()
            )
            pos.append(0)
        elif "n" in self.static_mods:
            mass.append(self.static_mods["n"])
            pos.append(0)

        # Extract the rest of the sequence
        residues = [res.group() for res in re.finditer("[A-Z](\[[^A-z]+\])?", seq)]
        for ind, res in enumerate(residues, 1):
            extracted_mod = re.search("(?<=\[)[^A-z]+(?=\])", res)
            if extracted_mod is not None:
                mass.append(STD_AA_MASS[res[0]] + float(extracted_mod.group()))
                pos.append(ind)
            elif res[0] in self.static_mods:
                mass.append(STD_AA_MASS[res[0]] + self.static_mods[res[0]])
                pos.append(ind)

        return np.array(pos, dtype=np.int32), np.array(mass, dtype=np.float32)


class IdentificationParser:
    """
    Parser for modification information coming from PSM file formats

    This class is designed to provide ease of access to PSMs from popular file
    formats. By providing the file of interest, the user can receive PSMs in
    either list form or dictionary form depending on need. PSM entries from
    individual file format are normalized to a single output format, which has
    the information necessary for running pyAscore. Sometimes modifications are
    not formatted in a straight forward way and there is usually a good chance
    that their masses will be trunctated. This class will attempt to normalize
    modifications so that all can be interpreted similarly.

    Parameters
    ----------
    id_file_name : str
        Path to file containing PSMs
    spec_file_format : str
        The format of the PSM containing file type. One of mzIdentML, pepXML,
        percolatorTXT, or mokapotTXT
    mass_corrector : MassCorrector
        Corrector class to normalize peptide modifications
    score_string : str
        String for score to extract from PSMs
    score_threshold : float
        Threshold to filter PSMs based on score
    score_lower_better : bool
        Whether a lower score is better than a higher score
    score_func : callable
        Transformation for scores
    static_mods : dict
        Dictionary of static mods to be used by some extractors.
    spec_file_name : str
        Currently not used
    """
    def __init__(self,
                 id_file_name, 
                 id_file_format,
                 mass_corrector=MassCorrector(),
                 score_string=None,
                 score_threshold=None,
                 score_lower_better=True,
                 score_func=None,
                 static_mods={"C": 57.021464},
                 spec_file_name=None):

        # Define a bridge to ID file handling classes
        if id_file_format == "mzIdentML":
            self._reader = MzIdentML(id_file_name, queue_size=32767)
            self._extractor = MzIdentMLExtractor(score_string)
        elif id_file_format == "pepXML":
            self._reader = PepXML(id_file_name, queue_size=32767)
            self._extractor = PepXMLExtractor(score_string)
        elif id_file_format == "percolatorTXT":
            self._reader = PercolatorTXT(id_file_name)
            self._extractor = PercolatorTXTExtractor(score_string, static_mods)
        elif id_file_format == "mokapotTXT":
            self._reader = MokapotTXT(id_file_name)
            self._extractor = MokapotTXTExtractor(score_string, static_mods)
        else:
            raise ValueError("{} not supported at this time."
                             " Must be on of: mzIdentML, pepXML,"
                             " percolatorTXT, or mokapotTXT".format(id_file_format))
        
        # Define containers for modification information
        # TODO: Allow for absolute modification masses
        self._match_records = []
        self.mass_corrector = mass_corrector

        # Define scoring information
        self.score_threshold = score_threshold
        self.score_lower_better = score_lower_better
        self.score_func = score_func

        # If tracking the spectra file name is desired, it is stored here
        # Right now this is here to allow building luciPhor input files 
        self.spec_file_name = spec_file_name


    def _get_match_records(self):
        """Internal function to build list of match records sorted by scan number

        If self._match_records is empty, this function fills the list, otherwise it does nothing.
        These records are unprocessed and will probably never be what the an end user wants, so 
        this function acts as a convenience to ensure that extraction happens once and only once.        
        """
        if not self._match_records:
            self._match_records = [e for e in self._reader.map(self._extractor.extract) if len(e['peptides']) > 0]
            self._match_records = sorted(self._match_records, key=lambda e: e["scans"][0])

    def _passes_scoring(self, score):
        """Check to see if PSM passes score threshold
        
        When no score_threshold is defined, this will always return True.

        Returns
        -------
        bool
            Whether a PSM should be retained
        """
        if self.score_threshold is None:
            return True

        elif score is None:
            return False

        else:
            return (score - self.score_threshold) * (-1 ** self.score_lower_better) > 0

    def _generate_hits(self):
        """Expands PSM entries and passes masses to be corrected
        
        Returns
        -------
        dict
            PSM object with normalized schema
        """
        self._get_match_records()
        for match in self._match_records:
            for ind in range(len(match['peptides'])):
                
                 
                mod_positions, mod_masses = self.mass_corrector.correct_multiple(match['peptides'][ind],
                                                                                 match['mod_positions'][ind],
                                                                                 match['mod_masses'][ind])

                score = match['scores'][ind]
                if self.score_func is not None and score is not None:
                    score = self.score_func(score)

                if not self._passes_scoring(score):
                    continue

                yield {"scan": match["scans"][ind],
                       "charge_state": match["charge_states"][ind],
                       "score": score,
                       "peptide": match['peptides'][ind],
                       "mod_positions": mod_positions,
                       "mod_masses": mod_masses}


    def to_list(self):
        """Return modified hits as list of dicts
        
        Returns
        -------
        list
            List of PSMs from file sorted by scan number
        """
        return [hit for hit in self._generate_hits()]

    def to_dict(self):
        """Return modified hits as dict

        Returns
        -------
        dict
            Dict of PSMs from file with schema: {scan number : spectra}
        """
        return {hit.pop("scan") : hit for hit in self._generate_hits()}
