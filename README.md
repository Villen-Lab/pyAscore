<img src="https://raw.githubusercontent.com/AnthonyOfSeattle/pyAscore/main/static/logo.png" width="300" title="pyAscore Logo">

## Intro
The **pyAscore** package provides a blazingly fast implementation of the Ascore algorithm for localizing peptide post-translational modifications from mass spectrometry data.
In order to provide efficient scoring, pyAscore implements dynamic programming over a custom modified peptide fragment tree and caches scoring calculations whenever possible.
This allows the algorithm to tackle both high and low resolution MS/MS spectra, as well as peptides of any length or number of modified amino acids.
The pyAscore package was also built without any assumptions on modification mass, and thus can be used to localize any feasible post-translational modification.
All algorithm components are implemented in C++, wrapped with Cython, and acessible by Python API or command line interface depending on pipeline needs.
