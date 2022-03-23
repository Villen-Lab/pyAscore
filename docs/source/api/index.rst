Python API
----------

The components of pyAscore can be imported directly into python in order to build
custom scripts for scoring the localization of algorithms. This allows advanced
logic to be performed such as spectral preprocessing or combining scores for more
confident localization. There are three main components that are described in the
subsequent pages:

   * **File reading:**
     Mass spectrometry data is usually disseminated in XML type files which can be
     read with the pyteomics python package. pyAscore provides some conveniece
     classes which build on this functionality and provide scan and PSM data in a
     convenient format for passing to scoring.

   * **Localization scoring:**
     The main scoring functionality of pyAscore is available as a class so that
     users have the option to tailor scoring to individual scans.

   * **Looping over modified peptides:** 
     Internally, pyAscore relies on an iterator which generates the theoretical
     fragments for every permutation of modified sites. We provide access to this
     iterator to fascilitate the building of new localization algorithms in python.

Example script
##############

The following is an example script which scores localization for PSMs in a single
mass spectrometry run. It runs almost exactly the same pipeline as the command
line version of pyAscore.

.. code-block:: python

   # 1) Read in PSMs
   psm_file = "psms.pep.xml"
   id_parser = IdentificationParser(psm_file, "pepXML")
   psm_objects = id_parser.to_list()
   
   # 2) Read in Spectra
   spectra_file = "spectra.mzML"
   spectra_parser = SpectraParser(spectra_file, "mzML")
   spectra_objects = spectra_parser.to_dict()

   # 3) Initialize pyAscore object
   mod_mass = 79.966331
   ascore = PyAscore(bin_size=100., n_top=10,
                     mod_group="STY",
                     mod_mass=mod_mass,
                     mz_error=.05)
   
   # 4) Score PSMs
   pyascore_results = []
   for psm in psm_objects:
       # 4.1) Check for modification of interest
       mod_select = np.isclose(psm["mod_masses"], mod_mass)
       nmods = np.sum(mod_select)
       
       if nmods >= 1:
           # 4.2) Grab spectrum
           spectrum = spectra_objects[psm["scan"]]

           # 4.3) Gather other modifications into aux mods
           aux_mod_pos = psm["mod_positions"][~mod_select].astype(np.uint32)
           aux_mod_masses = psm["mod_masses"][~mod_select].astype(np.float32)
           
           # 4.4) Run scoring algorithm
           ascore.score(mz_arr = spectrum["mz_values"],
                        int_arr = spectrum["intensity_values"],
                        peptide = psm["peptide"],
                        n_of_mod = np.sum(mod_select),
                        max_fragment_charge = psm["charge_state"] - 1,
                        aux_mod_pos = aux_mod_pos,
                        aux_mod_mass = aux_mod_masses)
           
           # 4.5) Place scores into an object to use later
           pyascore_results.append({"scan" : psm["scan"],
                                    "localized_peptide" : ascore.best_sequence,
                                    "pepscore" : ascore.best_score,
                                    "ascores" : ";".join([str(s) for s in ascore.ascores])})
   
   # 5) Make a dataframe of scores and write to a file                              
   pyascore_results = pd.DataFrame.from_records(pyascore_results)
   pyascore_results.to_csv("pyascore_results.tsv", sep="\t")

.. toctree::
   :maxdepth: 1
   :hidden:
   :titlesonly:

   Overview <self>
   parsing.rst
   scoring.rst
   modified_peptides.rst
