import re
import numpy as np
from numpy import isclose
from pyteomics.pepxml import PepXML
from pyteomics import mass

# Constants
STD_AA_MASS = mass.std_aa_mass

# In case a set of custom mods is not supplied
COMMON_MODS = {"n": 42.010565,
               "M": 15.9949,
               "S": 79.966331,
               "T": 79.966331,
               "Y": 79.966331,
               "C": 57.021464}

class MassCorrector:
    def __init__(self, mod_mass_dict=COMMON_MODS, aa_mass_dict=STD_AA_MASS, mz_tol=1.5, n_mod_ind=0):
        self.mod_mass_dict = mod_mass_dict
        self.aa_mass_dict = aa_mass_dict
        self.mz_tol = mz_tol
        self.n_mod_ind = n_mod_ind

    def correct(self, res, pos, mass):
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
            raise ValueError("Unrecognized mod on {} at position {} with mass: {}".format(res, pos, mass))

    def correct_multiple(self, peptide, positions, masses):
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


class IDExtractor:
    """
    Parser class for Pyteomics identifications dictionaries

    This class defines the core public functionality of identification extraction.
    Users will pass single entry dictionaries via the extract method, which will
    then call parsing methods specific to each file type, before sanitizing final
    results relavant to modifications, and passing these back in a common form.

    Attributes:
        score_string: Specific string describing the score to extract from a record.
        entry: Last Pyteomics entry parsed.
    """
    def __init__(self, score_string=None):
        """Initializes storage attributes for use by private methods."""
        self._match_list = None
        self._match = None
        self.score_string = score_string
        self.entry = None

    def _get_score(self):
        """Returns values requested via score_string if present"""
        try:
            return self._match[self.score_string]
        except KeyError:
            return None
    
    def _initialize_results(self):
        """
        Initializes the results dictionary with the correct sized list.
        
        Note: 
            I feel any speed increase from initializing the lists is low and makes
            the code less readable later.
        """
        nmatches = len(self._match_list)
        fields = ["scans", "scores", "charge_states",
                  "peptides", "mod_positions", "mod_masses"]
        self.results = {f: [None] * nmatches for f in fields}

    def extract(self, entry):
        """
        Sanatizes a single entry by calling filetype specific extraction methods

        Arguments:
            entry: A single Pyteomics dictionary of an identification file entry.
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

class PepXMLExtractor(IDExtractor):
    def _get_matches(self):
        try:
            return self.entry["search_hit"]
        except KeyError:
            return []

    def _get_scans(self):
        try:
            return int(self.entry["start_scan"])
        except KeyError:
            return -1

    def _get_charge(self):
        try:
            return int(self.entry["assumed_charge"])
        except KeyError:
            return 0

    def _get_score(self):
        try:
            return float(self._match["search_score"][self.score_string])
        except KeyError:
            return None

    def _get_peptide(self):
        try:
            return self._match["peptide"]
        except KeyError:
            return ""

    def _get_mod_info(self):
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


class IdentificationParser:
    """
    Parser for modification information coming from .pep.xmls.

    Sometimes modifications are not formatted in a straight forward way. There
    is usually a good chance that their masses will be trunctated or formatted
    differently than would be expected. This class is to make it easy to get
    sequences, modifications, and their masses easily from IDExtractor child
    classes.

    Attributes:
        
    """
    def __init__(self,
                 id_file_name, 
                 id_file_format,
                 mass_corrector=MassCorrector(), 
                 score_string=None,
                 score_threshold=None,
                 score_lower_better=True,
                 score_func=None,
                 spec_file_name=None):

        # Define a bridge to ID file handling classes
        if id_file_format == "pepXML":
            self._reader = PepXML(id_file_name)
            self._extractor = PepXMLExtractor(score_string)
        else:
            raise ValueError("{} not supported at this time.".format(id_file_format))
        
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
        """
        Internal function to build list of match records sorted by scan number

        If self._match_records is empty, this function fills the list, otherwise it does nothing.
        These records are unprocessed and will probably never be what the an end user wants, so 
        this function acts as a convenience to ensure that extraction happens once and only once.        
        """
        if not self._match_records:
            self._match_records = [e for e in self._reader.map(self._extractor.extract) if len(e['peptides']) > 0]
            self._match_records = sorted(self._match_records, key=lambda e: e["scans"][0])

    def _passes_scoring(self, score):
        if self.score_threshold is None:
            return True

        elif score is None:
            return False

        else:
            return (score - self.score_threshold) * (-1 ** self.score_lower_better) > 0

    def _generate_hits(self):
        """
        Expands PSM entries and passes masses to be corrected
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
        """
        Return modified hits as list of dicts with keys: (scan, charge_state, score, peptide, mod_positions, mod_masses)
        """
        return [hit for hit in self._generate_hits()]

    def to_dict(self):
        """
        Return modified hits as dict with scans as keys and dict values with keys: (charge_state, score, peptide, mod_positions, mod_masses)
        """
        return {hit.pop("scan") : hit for hit in self._generate_hits()}

#    def _format_mods(self, peptide, positions, masses):
#        mod_strings = []
#        for pos, mass in zip(positions, masses):
#            res = 'n' if pos == 0 else peptide[pos - 1]
#            res, pos, mass = self.mass_corrector.correct(res, pos, mass)
#            for p, m in zip(pos, mass):
#                mod_strings.append("{}={:.6f}".format(p if p == -100 else p-1, m))
#        return ",".join(mod_strings)
#
#    def _generate_tsv_entries(self):
#        self._get_match_records()
#        for match in self._match_records:
#            for ind in range(len(match['peptides'])):
#                mod_string = self._format_mods(match['peptides'][ind],
#                                               match['mod_positions'][ind],
#                                               match['mod_masses'][ind])
#                score = match['scores'][ind]
#                if self.score_func is not None and score is not None:
#                    score = self.score_func(score)
#
#                if len(mod_string) == 0 or not self._passes_scoring(score):
#                    continue
#
#                fields = [self.spec_file_name,
#                          match['scans'][ind],
#                          match['charge_states'][ind],
#                          score,
#                          match['peptides'][ind],
#                          mod_string]
#                yield "\t".join([str(e) for e in fields])
#                                     
#
#    def to_tsv(self, filename, sep="\t"):
#        with open(filename, "w") as dest:
#            dest.write("\t".join(["srcFile",
#                                  "scanNum",
#                                  "charge",
#                                  "PSMscore",
#                                  "peptide",
#                                  "modSites"]) + "\n")
#            for entry in self._generate_tsv_entries():
#                dest.write(entry + "\n")
#
