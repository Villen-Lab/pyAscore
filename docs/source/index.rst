Getting Started
---------------

**pyAscore** provides a fast and extensible version of the Ascore algorithm
for localizing post translational modifications (PTMs) from mass spectrometry
data. The package can be used from the command line, and reads many of the 
most popular data formats, or it can be incorporated into novel workflows
by using components directly in python.

Introduction
------------

Database search software tools such as Comet perform well at determining the
peptide sequence that generated a spectra and whether or not that peptide
contains a modification. Generally, they will place that modification on
one of the pre-specified amino acids, i.e. a phosphorylation will go on
a one of the S, T, or Y residues in a sequence, but they do not directly
tell you how confident they are in that placement.

In order to place a modification correctly on a peptide, you must observe
peaks which can only exist if the modification occurred at the specified 
site, i.e. site determining peaks. The Ascore algorithm was one of 
the first tools to explicitly score the confidence of PTM localization 
for a PSM. [1]_ It provides both the best localization for a PTM of 
interest on a peptide sequence as well as a probabilistic score for
how much better that localization is than the next best localization.

The pyAscore package provides a blazingly fast implementation of the Ascore
algorithm that can either be used from the command line or directly from
a Python script. The package's Cython backend implements dynamic programming
over a custom modified peptide fragment tree, and caches scoring calculations
whenever possible. This allows the algorithm to tackle both high and low resolution
MS/MS spectra, as well as peptides of any length or number of modified amino acids.
Furthermore, by specifying modification mass any PTM can be feasibly localized.

Please read bellow to see how you can install an use pyAscore, and checkout
our documentation to learn how you can apply our package to your own analyses.

Installation
------------

Before you can install and use pyAscore, you'll need to have Python 3.6+ and g++ 7+
installed. You can check versions for both with the following code:

.. code-block:: bash

   $ python3 --version
   $ g++ --version

pyAscore also depends on several Python packages:

- `numpy <https://numpy.org/>`_
- `scipy <https://scipy.org/>`_
- `pandas <https://pandas.pydata.org/>`_
- `cython <https://cython.org/>`_
- `lxml <https://lxml.de/>`_
- `pyteomics <https://pyteomics.readthedocs.io>`_
- `tqdm <https://github.com/tqdm/tqdm>`_

If you just want to use the pyAscore package, and don't want to contribute,
you can get the most up to date version using pip.

.. code-block:: bash

   $ pip install pyascore

If you would like to contribute, first fork the main repository, and then follow
the following steps to compile and test.

.. code-block:: bash

   $ git clone https://github.com/[USERNAME]/pyAscore.git
   $ cd pyAscore
   $ python setup.py build_ext --inplace
   $ python -m unittest 

Basic Usage
-----------

The main inputs to pyAscore are spectra from a masss spectrometry run and
PSMs from a database search of the spectra. For a full list of accepted formats
please see our command line interface page.

Run **pyAscore** from the command line
######################################

Once installed, the pyAscore package can be used straight from the command line.
By default the package will attempt to analyze the localization of phosphorylation,
but the modification of interest can be specified by the mass and the residues that
the modification can occupy. Analyses should also be tailored to the instrument, i.e.,
we recommend an `mz_error` of 0.05 for high resolution data and 0.5 for low resolution
data. A full list of parameters is available by running with the `-h` flag, or by
going to our command line interface page.

.. code-block:: bash

   $ pyascore --residues STY \
   >          --mod_mass 79.9663 \
   >          --mz_error .05 \
   >          spectra_file.mzML \
   >          psm_file.pep.xml \
   >          output_file.tsv


.. [1] Beausoleil, S. A., *et al.* "A probability-based approach for 
       high-throughput protein phosphorylation analysis and site localization.
       Nat. Biotechnol. 24, 1285–1292 (2006)


.. toctree::
   :hidden:
   :maxdepth: 1
   :caption: pyAscore
   :titlesonly:

   self
   cli.rst
   api/index.rst
