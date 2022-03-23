Scoring localization
--------------------

Here we will break down the scoring script from the first API page to give
an overview of the individual scoring components. The main scoring functionality
is performed on a PSM by PSM basis and can be accessed through the **pyAscore**
class. The format that I find most helpful in scripts is to read in the spectra
as a dictionary of dictionaries and the PSMs as a list so that you can loop through
the PSMs and then just look up the spectra you need.

.. code-block:: python

   psm_file = "psms.pep.xml"
   id_parser = IdentificationParser(psm_file, "pepXML")
   psm_objects = id_parser.to_list()

   spectra_file = "spectra.mzML"
   spectra_parser = SpectraParser(spectra_file, "mzML")
   spectra_objects = spectra_parser.to_dict()

The **pyAscore** class will score the modification of interest and has parameters
that help you tailor scoring to your individual instrument parameters. One thing
that is good to note is that the scoring is done with exact calculations that get
a lot of speed enhancement from caching of previous calculations. This means that
it is best to initialize ascore objects before looping through the PSMs, instead
of for every PSM.

.. code-block:: python

   mod_mass = 79.966331
   ascore = PyAscore(mod_group="STY",
                     mod_mass=mod_mass,
                     mz_error=.05,
                     fragment_types="by") 

The PSM reading functionality does not know which modification is of interest to
the user, and so during looping the modifications must be partitioned into variable
and static mods. Then, the fragments, intensities, unmodified peptide sequence, and
modification information for a PSM can be passed to the score function. At this point,
the maximum fragment charge state can also be chosen. For max speed, we can only score
+1 peaks, but for max accuracy we would recommend up to the precursor charge - 1.

.. code-block:: python

   pyascore_results = []
   for psm in psm_objects:
       mod_select = np.isclose(psm["mod_masses"], mod_mass)
       nmods = np.sum(mod_select)

       if nmods >= 1:
           spectrum = spectra_objects[psm["scan"]]
           
           # Partition modifications
           aux_mod_pos = psm["mod_positions"][~mod_select].astype(np.uint32)
           aux_mod_masses = psm["mod_masses"][~mod_select].astype(np.float32)

           # Score PSMs
           ascore.score(mz_arr = spectrum["mz_values"],
                        int_arr = spectrum["intensity_values"],
                        peptide = psm["peptide"],
                        n_of_mod = np.sum(mod_select),
                        max_fragment_charge = psm["charge_state"] - 1,
                        aux_mod_pos = aux_mod_pos,
                        aux_mod_mass = aux_mod_masses)

           # Store scores
           pyascore_results.append({"scan" : psm["scan"],
                                    "localized_peptide" : ascore.best_sequence,
                                    "pepscore" : ascore.best_score,
                                    "ascores" : ";".join([str(s) for s in ascore.ascores])})


Class Reference
###############

.. autoclass:: pyascore.PyAscore
   :members:

