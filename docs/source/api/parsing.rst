Reading data files
------------------

pyAscore relies on the functionality of pyteomics to read XML files
of spectra and PSMs, with tsv file reading being provided by pandas.
While this system is a powerful way to extract information from files,
pyAscore also provides some convenience functions to extract just the
relevant information for scoring localizations from files in a standard
format. This functionality can be accessed through the **SpectraParser**
and **IdentificationParser** classes.

Reading spectra
###############

Spectra file reading is fairly straight forward, and can be achieved by
specifying the file name and file type. An optional argument for MSn 
level is also provided, but the default MSn level of 2 should be correct
for most purposes. Spectra are read in parallel if multiple cores are
available and can be transformed to a list of dictionaries or a dictionary
of dictionaries.

.. code-block:: python

   spectra_file = "spectra.mzML"
   spectra_parser = SpectraParser(spectra_file, "mzML")

   # To a list of dictionaries:
   spectra_list = spectra_parser.to_list()

   # To a dictionary of dictionaries:
   spectra_dict = spectra_parser.to_dict()

Reading PSMs
############

Reading PSM files happens in much the same way as spectra by specifying a
file name and the file type.

.. code-block:: python

   psm_file = "psms.pep.xml"
   psm_parser = IdentificationParser(psm_file, "pepXML")

   # To a list of dictionaries:
   psm_list = psm_parser.to_list()

   # To a dictionary of dictionaries:
   psm_dict = psm_parser.to_dict()

Likely, the most useful extra feature of this module is the mass correction
that allows you to make sure that the correct modification mass is associated
with residues when they are read in. This is important since many search
engines and other programs will truncate the mass. This functionality is
supplied by the **MassCorrector** class. It comes with several modifications
built in, but if you happen to have one that isn't recognized, then it would
be good to use the following code.

.. code-block:: python

   modifications = {"n": 42.010565, # N-term acetylation
                    "M": 15.9949,   # Methionine oxidation
                    "S": 79.966331, # Serine Phoshorylation
                    "T": 79.966331, # Threonine Phosphorylation
                    "Y": 79.966331, # Tyrosine Phosphorylation
                    "C": 57.021464} # Cysteine Carbamidomethylation
   mass_corrector = MassCorrector(modifications, mz_tol=1.5)

   psm_file = "psms.pep.xml"
   psm_parser = IdentificationParser(psm_file, 
                                     "pepXML",
                                     mass_corrector)

Class Reference
###############

.. autoclass:: pyascore.SpectraParser
   :members:

.. autoclass:: pyascore.IdentificationParser
   :members:

.. autoclass:: pyascore.MassCorrector

