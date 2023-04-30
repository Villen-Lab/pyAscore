<img src="https://raw.githubusercontent.com/AnthonyOfSeattle/pyAscore/main/static/logo.png" width="300" title="pyAscore Logo">

[![Linux Build](https://github.com/Villen-Lab/pyAscore/actions/workflows/linux-build.yml/badge.svg)](https://github.com/Villen-Lab/pyAscore/actions/workflows/linux-build.yml)
[![Mac Build](https://github.com/Villen-Lab/pyAscore/actions/workflows/mac-build.yml/badge.svg)](https://github.com/Villen-Lab/pyAscore/actions/workflows/mac-build.yml)
[![Windows Build](https://github.com/Villen-Lab/pyAscore/actions/workflows/windows-build.yml/badge.svg)](https://github.com/Villen-Lab/pyAscore/actions/workflows/windows-build.yml)
[![Documentation Status](https://readthedocs.org/projects/pyascore/badge/?version=latest)](https://pyascore.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/pyascore?color=green)](https://pypi.org/project/pyascore/)

## Intro

The **pyAscore** package provides a blazingly fast implementation of the Ascore algorithm for localizing peptide post-translational modifications from mass spectrometry data.
In order to provide efficient scoring, pyAscore implements dynamic programming over a custom modified peptide fragment tree and caches scoring calculations whenever possible.
This allows the algorithm to tackle both high and low resolution MS/MS spectra, as well as peptides of any length or number of modified amino acids.
The pyAscore package was also built without any assumptions on modification mass, and thus can be used to localize any feasible post-translational modification.
All algorithm components are implemented in C++, wrapped with Cython, and acessible by Python API or command line interface depending on pipeline needs.

For more information, check out our
[documentation](https://pyascore.readthedocs.io).

## Getting Started

### Install from PyPI

If you just want to use the pyAscore package, and don't want to contribute, you can get the most up to date version using pip.

```
pip install pyascore
```

### Installing from a local clone

If you would like to contribute, first fork the main repository, and then follow the following steps to compile and test.

```
git clone https://github.com/[USERNAME]/pyAscore.git
cd pyAscore
python setup.py build_ext --inplace
python -m unittest 
```
### Usage

The pyAscore package can be used straight from the command line as a module. 
A full list of parameters is available by running with the `-h` flag.

```
$ pyascore -h

usage: pyAscore [-h] [--match_save] [--residues RESIDUES]
                [--mod_mass MOD_MASS] [--mz_error MZ_ERROR]
                [--mod_correction_tol MOD_CORRECTION_TOL]
                [--zero_based ZERO_BASED]
                [--neutral_loss_groups NEUTRAL_LOSS_GROUPS]
                [--neutral_loss_masses NEUTRAL_LOSS_MASSES]
                [--static_mod_groups STATIC_MOD_GROUPS]
                [--static_mod_masses STATIC_MOD_MASSES]
                [--fragment_types FRAGMENT_TYPES]
                [--max_fragment_charge MAX_FRAGMENT_CHARGE]
                [--hit_depth HIT_DEPTH] [--parameter_file PARAMETER_FILE]
                [--spec_file_type SPEC_FILE_TYPE]
                [--ident_file_type IDENT_FILE_TYPE]
                spec_file ident_file out_file

The pyAscore module provides PTM localization analysis using a custom
implementation of the Ascore algorithm. It employees pyteomics for efficient
reading of spectra in mzML format and identifications in pepXML format. All
scoring has been implemented in custom c++ code which is exposed to python via
cython wrappers. Any PTM which be defined with a canonical amino acid and mass
shift can be analyzed.

positional arguments:
  spec_file             MS Spectra file.
  ident_file            Results of database search.
  out_file              Destination for Ascores.

optional arguments:
  -h, --help            show this help message and exit
  --match_save
  --residues RESIDUES   Residues which can be modified.
  --mod_mass MOD_MASS   Modification mass to match to identifications. This is
                        often rounded by search engines so this argument
                        should be considered the most accurate mass.
  --mz_error MZ_ERROR   Tolerance in mz for deciding whether a spectral peak
                        matches to a theoretical peak.
  --mod_correction_tol MOD_CORRECTION_TOL
                        MZ tolerance for deciding whether a reported
                        modification matches internal or user specified
                        modifications. A wide tolerance can help overcome
                        rounding. If more precission is needed, make sure to
                        set this parameter and that your search engine
                        provides for it.
  --zero_based ZERO_BASED
                        Mod positions are by default assumed to be 1 based.
  --neutral_loss_groups NEUTRAL_LOSS_GROUPS
                        Comma separated clusters of amino acids which are
                        expected to have a neutral loss. To specify that the
                        modified versions of the amino acids should have the
                        neutral loss, use lower case letters. Example: 'st' vs
                        'ST'.
  --neutral_loss_masses NEUTRAL_LOSS_MASSES
                        Comma separated neutral loss masses for each of the
                        neutral_loss_groups. Should have one mass per group.
                        Positive masses indicate a loss, e.g. '18.0153' for
                        water loss, while negative masses can be used to
                        indicate a gain.
  --static_mod_groups STATIC_MOD_GROUPS
                        Comma separated clusters of amino acids which will be
                        read in with a constant modification.
  --static_mod_masses STATIC_MOD_MASSES
                        Comma separated masses for each of the
                        static_mod_groups.
  --fragment_types FRAGMENT_TYPES
                        Fragment ion types to score. Supported: bcyzZ. The
                        special character Z indicates a z+H fragment.
  --max_fragment_charge MAX_FRAGMENT_CHARGE
                        Max fragment charge to use for calculating theoretical
                        peaks. Internally, the max fragment charge will not be
                        allowed to be greater than the PSM charge - 1.
                        However, if a more stringent limit needs to be set,
                        this argument can be used.
  --hit_depth HIT_DEPTH
                        Number of PSMS to take from each scan. Set to negative
                        to always analyze all.
  --parameter_file PARAMETER_FILE
                        A file containing parameters. Example: 'residues =
                        STY'.
  --spec_file_type SPEC_FILE_TYPE
                        The type of file supplied for spectra. One of mzML or
                        mzXML. Default: mzML.
  --ident_file_type IDENT_FILE_TYPE
                        The type of file supplied for identifications. One of
                        pepXML, mzIdentML, percolatorTXT, or mokapotTXT.
                        Default: pepXML.
```
