import re
import numpy as np
from pyteomics.mzml import MzML
from pyteomics.mzxml import MzXML


class SpectraExtractor:
    def __init__(self):
        self.scan = None
        self.results = None

    def extract(self, scan):
        self.scan = scan
        
        self.results = {}
        self.results["scan"] = self._get_scan()
        self.results["ms_level"] = self._get_ms_level()
        self.results["precursor_mz"],  self.results["precursor_charge"] = self._get_precursor()
        self.results["mz_values"], self.results["intensity_values"] = self._get_spectra()

        return self.results


class MzMLExtractor(SpectraExtractor):
    def _get_scan(self):
        try:
            scan_string = re.search("(?<=scan=)[0-9]+", self.scan["id"]).group()
            return int(scan_string)
        except KeyError:
            return -1

    def _get_ms_level(self):
        try:
            return self.scan["ms level"]
        except KeyError:
            return 0

    def _get_precursor(self):
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
        try:
            return self.scan["m/z array"], self.scan["intensity array"]
        except:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)


class MzXMLExtractor(SpectraExtractor):
    def _get_scan(self):
        try:
            return int(self.scan["num"])
        except KeyError:
            return -1

    def _get_ms_level(self):
        try:
            return int(self.scan["msLevel"])
        except KeyError:
            return 0

    def _get_precursor(self):
        try:
            precursor_list = self.scan["precursorMz"]
            nprecursors = len(precursor_list)
            if nprecursors > 1:
                 raise ValueError("Multiple precursors not supported at this time")

            return precursor_list[0]["precursorMz"], precursor_list[0]["precursorCharge"]

        except KeyError:
            return None, None

    def _get_spectra(self):
        try:
            return self.scan["m/z array"], self.scan["intensity array"]
        except:
            return np.array([], dtype=np.float64), np.array([], dtype=np.float64)


class SpectraParser:
    """
    Parser for spectra information coming from mzMLs or mzXMLs.
    """
    def __init__(self,
                 spec_file_name,
                 spec_file_format,
                 ms_level=2,
                 custom_filter=None):

        # Define a bridge to spectra file handling classes
        if spec_file_format == "mzML":
            self._reader = MzML(spec_file_name)
            self._extractor = MzMLExtractor()
        elif spec_file_format == "mzXML":
            self._reader = MzXML(spec_file_name)
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
        passes_ms_level = True
        if self.ms_level:
            passes_ms_level = entry["ms_level"] == self.ms_level

        passes_custom_filter = True
        if self.custom_filter is not None:
            passes_custom_filter = self.custom_filter(entry)

        return passes_ms_level and passes_custom_filter

    def _get_spectra(self):
        """
        Internal function to build list of scans sorted by scan number       
        """
        if not self._spectra:
            self._spectra = [s for s in self._reader.map(self._extractor.extract) if self._passes_filtering(s)]
            self._spectra = sorted(self._spectra, key=lambda s: s["scan"])

    def to_list(self):
        self._get_spectra()
        return self._spectra

    def to_dict(self):
        self._get_spectra()
        return {spec.pop("scan") : spec for spec in self._spectra} 
