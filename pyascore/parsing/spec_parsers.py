import re
import numpy as np
from pyteomics.mzml import MzML
from pyteomics.mzxml import MzXML


class SpectraExtractor:
    """
    Abstract class providing a common callable to transform scans from pyteomics

    Different spectra file formats provide scan objects with specific schemas that
    are not necessarily cross compatible. To make it easy to take scans from any
    file format, classes that inherit from this class provide a series of methods
    to return specific attributes from scans. This allows the class to transform
    scans into a common schema containing only the information necessary for scoring.
    """
    def __init__(self):
        self.scan = None
        self.results = None

    def extract(self, scan):
        """Extract information from a pyteomics scan dictionary

        This function can be passed to the pyteomics map method to extract useful
        informations from scans. Derivitive classes call internal methods to
        retrieve specific information in file format specific queries.

        Parameters
        ----------
        scan : dict
            Scan dictionary passed by pyteomics

        Returns
        -------
        dict
            Scan information with common schema
        """
        self.scan = scan
        
        self.results = {}
        self.results["scan"] = self._get_scan()
        self.results["ms_level"] = self._get_ms_level()
        self.results["precursor_mz"],  self.results["precursor_charge"] = self._get_precursor()
        self.results["mz_values"], self.results["intensity_values"] = self._get_spectra()

        return self.results


class MzMLExtractor(SpectraExtractor):
    def _get_scan(self):
        """Extract scan number from scan

        Returns
        -------
        int
            Scan number
        """
        try:
            scan_string = re.search("(?<=scan=)[0-9]+", self.scan["id"]).group()
            return int(scan_string)
        except KeyError:
            return -1

    def _get_ms_level(self):
        """Extract MSn level from scan, normally 1 or 2

        Returns
        -------
        int
            MSn level
        """
        try:
            return self.scan["ms level"]
        except KeyError:
            return 0

    def _get_precursor(self):
        """Extract precursor m/z and charge state from scan

        Returns
        -------
        (float, int)
            Tuple with precursor m/z (0) and charge state (1)
        """
        try:
            precursor_list = self.scan["precursorList"]
            nprecursors = precursor_list["count"]
            if nprecursors > 1:
                 raise ValueError("Multiple precursors not supported at this time")
            
            ion_dict = precursor_list["precursor"][0]["selectedIonList"]["selectedIon"][0]
            return ion_dict["selected ion m/z"], ion_dict["charge state"]

        except KeyError:
            return None, None

    def _get_spectra(self):
        """Extract fragment m/z array and intensity array from scan

        Returns
        -------
        (np.ndarray, np.ndarray)
            Tuple with fragment m/z array (0) and intensity array (1)
        """
        try:
            mzs = self.scan["m/z array"].astype(np.float64)
            intensities = self.scan["intensity array"].astype(np.float64)
            return mzs, intensities
        except:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)


class MzXMLExtractor(SpectraExtractor):
    def _get_scan(self):
        """Extract scan number from scan

        Returns
        -------
        int
            Scan number
        """
        try:
            return int(self.scan["num"])
        except KeyError:
            return -1

    def _get_ms_level(self):
        """Extract MSn level from scan, normally 1 or 2

        Returns
        -------
        int
            MSn level
        """
        try:
            return int(self.scan["msLevel"])
        except KeyError:
            return 0

    def _get_precursor(self):
        """Extract precursor m/z and charge state from scan

        Returns
        -------
        (float, int)
            Tuple with precursor m/z (0) and charge state (1)
        """
        try:
            precursor_list = self.scan["precursorMz"]
            nprecursors = len(precursor_list)
            if nprecursors > 1:
                 raise ValueError("Multiple precursors not supported at this time")

            return precursor_list[0]["precursorMz"], precursor_list[0]["precursorCharge"]

        except KeyError:
            return None, None

    def _get_spectra(self):
        """Extract fragment m/z array and intensity array from scan

        Returns
        -------
        (np.ndarray, np.ndarray)
            Tuple with fragment m/z array (0) and intensity array (1)
        """
        try:
            mzs = self.scan["m/z array"].astype(np.float64)
            intensities = self.scan["intensity array"].astype(np.float64)
            return mzs, intensities
        except:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)


class SpectraParser:
    """
    Parser to read spectra from mzML and mzXML files

    This class is designed to provide ease of access to spectra from popular file formats.
    By providing the file of interest and the MSn level, the user can receive spectra in
    either list form or dictionary form depending on need. Spectral entries from individual
    file formats are normalized to a single output schema, which has the information
    necessary for running pyAscore. This output can be further filtered by supplying a
    custom filter, and future versions will make this option more powerful.

    Parameters
    ----------
    spec_file_name : str
        Path to spectral file
    spec_file_type : str
        The spectra file's type. One of mzML or mzXML
    ms_level : int
        MSn level to be returned to the user
    custom_filter : callable
        A callable which takes a spectral object and returns a boolean which states whether
        a spectra should be retained. The spectral objects passed to this parameter are
        currently the same as the ones returned to the user
    """
    def __init__(self,
                 spec_file_name,
                 spec_file_format,
                 ms_level=2,
                 custom_filter=None):

        # Define a bridge to spectra file handling classes
        if spec_file_format == "mzML":
            self._reader = MzML(spec_file_name, queue_size=32767)
            self._extractor = MzMLExtractor()
        elif spec_file_format == "mzXML":
            self._reader = MzXML(spec_file_name, queue_size=32767)
            self._extractor = MzXMLExtractor()
        else:
            raise ValueError("{} not supported at this time."
                             " Should be one of: mzML or mzXML".format(spec_file_format))


        # Store filtering criteria
        if ms_level >= 0:
            self.ms_level = ms_level
        else:
            raise ValueError("ms_level must be an integer greater than or equal to 0")

        self.custom_filter = None
        if custom_filter is not None:
            if callable(ms_level):
                self.custom_filter = custom_filter
            else:
                raise ValueError("custom_filter must be callable.")

        # Storage
        self._spectra = []

    def _passes_filtering(self, entry):
        """Internal function to decide whether a spectra should be retained
        
        Returns
        -------
        bool
            True if scan should be retained, else False
        """

        passes_ms_level = True
        if self.ms_level:
            passes_ms_level = entry["ms_level"] == self.ms_level

        passes_custom_filter = True
        if self.custom_filter is not None:
            passes_custom_filter = self.custom_filter(entry)

        return passes_ms_level and passes_custom_filter

    def _get_spectra(self):
        """Internal function to build list of scans sorted by scan number"""

        if not self._spectra:
            self._spectra = [s for s in self._reader.map(self._extractor.extract)
                             if self._passes_filtering(s)]
            self._spectra = sorted(self._spectra, key=lambda s: s["scan"])

    def to_list(self):
        """Return spectra from file in list form
        
        Returns
        -------
        list
            List of scans from file sorted by scan number
        """

        self._get_spectra()
        return self._spectra

    def to_dict(self):
        """Return spectra from file in dictionary form
        
        Returns
        -------
        dict
            Dict of scans from file with schema: {scan number : spectra}
        """

        self._get_spectra()
        return {spec.pop("scan") : spec for spec in self._spectra} 
