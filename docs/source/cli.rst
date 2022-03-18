Command Line Interface
######################

.. contents::
   :depth: 2
   :local:

Common workflows
----------------

Running with standard file formats
==================================

The minium information needed to run pyAscore from the command line
is a file containing spectra and a file containing PMSs from a database
search engine such as Comet. By default, pyAscore will accept spectra
in the mzML format and PSMs in pepXML format. The modification of interest
can be specified by providing the modifiable amino acids in there 1 letter
codes and the mass fo the modification, e.g. S (Serine), T (Threonine), and 
Y (Tyrosine) and 79.9663 Da. This is also the default modification.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          spectra_file.mzML \
   >          psm_file.pep.xml \
   >          output_file.tsv

Other file formats can also be provided. For spectra, users can also supply
mzXML files, and for PSMs, users can supply mzIdentML, percolatorTXT, and
mokapotTXT. The format of the input files must be specified with the
`--spec_file_type` and `--ident_file_type` arguments respectively.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          --spec_file_type mzXML \
   >          --ident_file_type mzIdentML \
   >          spectra_file.mzXML \
   >          psm_file.mzid \
   >          output_file.tsv

Inputing custom PSM data
========================

Often, a user may have used a search engine which doesn't output a standard
format but they still want to connect it to pyAscore. Or they may want to
manipulate data before handing it off to pyAscore and want to work with a
simpler format for connecting their pipeline than XML based data. This
situation is easily handled by the percolatorTXT input format since the
standard Percolator output can be reduced to a minimal tab delimited layout.

.. code-block:: bash

      scan  rt     sequence                      charge
   0  2082  304.8  S[79.9663]NNSNSNSGGK          2     
   1  2624  402.8  RARES[79.9663]DNEDAK          3     
   2  2625  402.9  SS[79.9663]NGNESNGAK          2     
   3  2655  405.4  n[42.010565]S[79.9663]DAGRK   2     

Once a user makes a tsv file with their PSMs, the data can be fed to
pyAscore with the following command.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          --ident_file_type percolatorTXT \
   >          spectra_file.mzML \
   >          psm_file.tsv \
   >          output_file.tsv

Tailoring to instrument parameters
==================================

pyAscore works by looking for fragment ions within the supplied spectra
which match a theoretical fragment pattern. Users should tailor the
tolerance of this search so that it matches the instrument resolution.
This can be done with the `--mz_error` option. For high resolution data,
we have been using a tolerance of 0.05 Da, and for low resolution data,
we have been using a tolerance of 0.5 Da.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          --mz_error 0.05 \
   >          spectra_file.mzML \
   >          psm_file.pep.xml \
   >          output_file.tsv

By default, pyAscore will score `b` and `y` ion fragment peaks, which
are the most abundant peaks for HCD and CID fragmentation data. If a
user wants to analyze ETD fragmentation data, it is recommended to score
`c` and `z+H` ions. This can be specified with the `--fragment_types` option.
The `Z` character is used to differentiate the `z+H` ion from the `z` ion.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          --mz_error 0.05 \
   >          --fragment_types cZ \
   >          spectra_file.mzML \
   >          psm_file.pep.xml \
   >          output_file.tsv

Specifying neutral losses
=========================

A user has the option to use neutral loss peaks in their scoring procedure,
and these can be different for modified and unmodified residues. A user can
supply a comma sepparated list of amino acid groups, uppercase for unmodified
residues and lowercase for modified, and a comma sepparated list of neutral
loss masses. If, for example, a user wants to use the H3P04 neutral loss ions,
a loss of 97.976896, on modified Ser, Thr, and Tyr residues, they could use
the following command.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          --mz_error 0.05 \
   >          --neutral_loss_groups sty \
   >          --neutral_loss_masses 97.976896 \
   >          spectra_file.mzML \
   >          psm_file.pep.xml \
   >          output_file.tsv

While a rare occurence, a user could theoretically specify a gain of mass
on any residue by passing negative masses to `--neutral_loss_masses`.

Output description
------------------

pyAscore outputs a single `.tsv` style file with one entry for every
PSM containing the modification of interest. Example output is given
below.

.. code-block:: bash

   Scan	LocalizedSequence	                         PepScore           Ascores	                  AltSites
   7546	ARS[80]VS[80]PPPK	                         inf	            inf;inf	                  ;
   7547	S[80]AS[80]SC[57]PNLLVPETWPHQVSASHAGRSKQP        6.8605122566223145 0.0;0.0	                  4;4
   7548	VGSLM[16]TSSSGTSLRTSST[80]	                 16.59139633178711  0.0	                          11,12,15,16,17
   7549	NDSLSSLDFDDDDVDLS[80]REK	                 2.4440932273864746 8.094615	                  3,5,6
   7552	ASAS[80]PSTSSTSSRPK	                         92.12223815917969  0.0	                          2
   7553	RLNHS[80]PPQSSSR	                         31.98526954650879  53.332687	                  9,10,11
   7555	M[16]HSGEKPY[80]EC[57]S[80]EC[57]GKIFS[80]M[16]K 8.363080978393555  0.0;0.0;15.29163	          3;3;3
   7557	RHS[80]HS[80]HS[80]PMSTR	                 66.47030639648438  44.625744;48.631317;39.542572 10,11;10,11;10,11

Description of columns:
   * **Scan:**
     This is the scan number from the supplied spectra file.
     It is usually taken from the scan header and so care should
     by taken that this matches expectations.
   * **LocalizedSequence:**
     This is the modified peptide with the PTM of interest placed
     in the best positions according the the PepScore. All outputed
     masses are rounded to their whole number representations.
   * **PepScore:**
     This score gives the total amount of evidence for the listed
     sequence being correct. It is based on the total number of
     matching theoretical ions to the ranked peaks within the
     supplied spectrum. A value of `inf` means that there
     is no ambiguity in the localized sequence.
   * **Ascores:**
     This semicolon separated list of scores gives the relative
     amount of evidence for the localization of the modification
     of interest vs the next best localization. It is based on
     the number of matching theoretical site determining ion peaks
     in the supplied spectrum. A value of `inf` means that there
     is no ambiguity in the site placement. There is one entry in
     the list per modification of interest on the peptide.
   * **Altsites:**
     This semicolon separated list of comma separated positions
     gives the next best locations for a modification. There is one
     list of alternative sites per modification of interest on the peptide.

All options
-----------

.. argparse::
   :module: pyascore.config
   :func: build_parser
   :prog: pyascore

