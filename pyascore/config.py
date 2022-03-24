import re
import sys
import argparse

def args_from_file(file_path):
    """Read file with parameters listed in param=value pairs"""
    
    with open(file_path, "r") as parameters:
        arg_list = []
        for line in parameters:
            line = line.split("#")[0].rstrip().lstrip()
            search_result = re.search("^([^\s]+)\s*=\s*([^\s]+)$", line)
            if search_result is not None:
                arg_list.append("--" + search_result.group(1))
                arg_list.append(search_result.group(2))

        return arg_list

def build_parser():
    """Build parser for command line arguments."""

    parser = argparse.ArgumentParser(
        prog="pyAscore", description="The pyAscore module provides PTM localization analysis"
                                     " using a custom implementation of the Ascore algorithm."
                                     " It employees pyteomics for efficient reading of spectra"
                                     " in mzML format and identifications in pepXML format."
                                     " All scoring has been implemented in custom c++ code which"
                                     " is exposed to python via cython wrappers."
                                     " Any PTM which be defined with a canonical amino acid and"
                                     " mass shift can be analyzed."
    )
    parser.add_argument("--match_save", action="store_true")
    parser.add_argument("--residues", type=str, default="STY",
                        help="Residues which can be modified.")
    parser.add_argument("--mod_mass", type=float, default=79.966331,
                        help="Modification mass to match to identifications."
                             " This is often rounded by search engines so this"
                             " argument should be considered the most accurate mass.")
    parser.add_argument("--mz_error", type=float, default=0.5,
                        help="Tolerance in mz for deciding whether a spectral peak"
                             " matches to a theoretical peak.")
    parser.add_argument("--mod_correction_tol", type=float, default=1.,
                        help="MZ tolerance for deciding whether a reported modification"
                             " matches internal or user specified modifications."
                             " A wide tolerance can help overcome rounding."
                             " If more precission is needed, make sure to"
                             " set this parameter and that your search"
                             " engine provides for it.")
    parser.add_argument("--zero_based", type=bool, default=False,
                        help="Mod positions are by default assumed to be 1 based.")
    parser.add_argument("--neutral_loss_groups", type=str, default="",
                        help="Comma separated clusters of amino acids"
                             " which are expected to have a neutral loss."
                             " To specify that the modified versions of the amino acids"
                             " should have the neutral loss, use lower case letters."
                             " Example: 'st' vs 'ST'.")
    parser.add_argument("--neutral_loss_masses", type=str, default="",
                        help="Comma separated neutral loss masses for each of the neutral_loss_groups."
                             " Should have one mass per group."
                             " Positive masses indicate a loss, e.g. '18.0153' for water loss,"
                             " while negative masses can be used to indicate a gain.")
    parser.add_argument("--static_mod_groups", type=str, default="C",
                        help="Comma separated clusters of amino acids"
                             " which will be read in with a constant modification.")
    parser.add_argument("--static_mod_masses", type=str, default="57.021464",
                        help="Comma separated masses for each of the static_mod_groups.")
    parser.add_argument("--fragment_types", type=str, default="by",
                        help="Fragment ion types to score. Supported: bcyzZ."
                             " The special character Z indicates a z+H fragment.")
    parser.add_argument("--max_fragment_charge", type=int, default=5,
                        help="Max fragment charge to use for calculating theoretical peaks."
                             " Internally, the max fragment charge will not be allowed to be"
                             " greater than the PSM charge - 1. However, if a more stringent"
                             " limit needs to be set, this argument can be used.")
    parser.add_argument("--hit_depth", type=int, default=1,
                        help="Number of PSMS to take from each scan."
                           " Set to negative to always analyze all.")
    parser.add_argument("--parameter_file", type=str, default="",
                        help="A file containing parameters. Example: 'residues = STY'.")
    parser.add_argument("--spec_file_type", type=str, default="mzML",
                        help="The type of file supplied for spectra."
                             " One of mzML or mzXML. Default: mzML.")
    parser.add_argument("--ident_file_type", type=str, default="pepXML",
                        help="The type of file supplied for identifications."
                             " One of pepXML, mzIdentML, percolatorTXT, or mokapotTXT. Default: pepXML.")
    parser.add_argument("spec_file", type=str,
                        help="MS Spectra file.")
    parser.add_argument("ident_file", type=str,
                        help="Results of database search.")
    parser.add_argument("out_file", type=str,
                        help="Destination for Ascores.")

    return parser
