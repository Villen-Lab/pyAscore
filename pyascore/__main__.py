import sys
import argparse
import re
from pyascore import *
from datetime import datetime
from itertools import groupby
import numpy as np
from pandas import DataFrame
from tqdm import tqdm
import pickle

def get_time_stamp():
    return datetime.now().strftime("%m/%d/%y %H:%M:%S")

def args_from_file(file_path):
    with open(file_path, "r") as parameters:
        arg_list = []
        for line in parameters:
            line = line.split("#")[0].rstrip().lstrip()
            search_result = re.search("^([^\s]+)\s*=\s*([^\s]+)$", line)
            if search_result is not None:
                arg_list.append("--" + search_result.group(1))
                arg_list.append(search_result.group(2))
        
        return arg_list

def build_spectra_parser(arg_ref):
    print("{} -- Reading spectra from: {}".format(get_time_stamp(), arg_ref.spec_file))
    spec_parser = SpectraParser(arg_ref.spec_file, "mzML")
    return spec_parser.to_dict()


def build_identification_parser(arg_ref):
    print("{} -- Reading identifications from: {}".format(get_time_stamp(), arg_ref.ident_file))
    id_parser = iter(
        sorted(IdentificationParser(arg_ref.ident_file, "pepXML").to_list(),
               key=lambda spec: spec["scan"])
    )
    return id_parser


def build_ascore(arg_ref):
    return PyAscore(bin_size=100., n_top=10,
                    mod_group=arg_ref.residues,
                    mod_mass=arg_ref.mod_mass,
                    mz_error=arg_ref.mz_error)

def process_mods(arg_ref, positions, masses):
    variable_mod_count = 0
    const_pos, const_masses = [], []
    for pos, mass in zip(positions, masses):
        if np.isclose(arg_ref.mod_mass, mass,
                      rtol=1e-6, atol=arg_ref.mod_correction_tol):
            variable_mod_count += 1
        else:
            shift = 1 if not arg_ref.zero_based else 0
            const_pos.append(pos - shift)
            const_masses.append(mass)

    return (np.array(const_pos, dtype=np.uint32),
            np.array(const_masses, dtype=np.float32),
            variable_mod_count)


def save_match(spectra, match):
    with open("dump_spectra.pkl", "wb") as src:
        pickle.dump([spectra], src)

    with open("dump_match.pkl", "wb") as src:
        pickle.dump([match], src)


def main():
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
                        help="Residues which can be modified")
    parser.add_argument("--mod_mass", type=float, default=79.966331,
                        help="Modification mass to match to identifications."
                             " This is often rounded by search engines so this"
                             " argument should be considered the most accurate mass")
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
                        help="Mod positions are by default assumed to be 1 based")
    parser.add_argument("--hit_depth", type=int, default=1,
                        help="Number of PSMS to take from each scan."
                           " Set to negative to always analyze all.")
    parser.add_argument("--parameter_file", type=str, default="",
                        help="A file containing parameters. e.x. param = val")
    parser.add_argument("spec_file", type=str,
                        help="MS Spectra file supplied as MzML")
    parser.add_argument("ident_file", type=str,
                        help="Comet hits supplied as pepXML")
    parser.add_argument("out_file", type=str,
                        help="Destination for Ascores")
    args = parser.parse_args()

    if args.parameter_file:
        args = parser.parse_args(args_from_file(args.parameter_file) + sys.argv[1:])

    print("{} -- Ascore Started".format(get_time_stamp()))

    spectra_map = build_spectra_parser(args)
    psms = build_identification_parser(args)

    print("{} -- Anlyzing PSMs".format(get_time_stamp()))
    ascore = build_ascore(args)
    scores = []
    for _, match_group in tqdm(groupby(psms, lambda x: x["scan"])):
        for match_ind, match in enumerate(match_group):
            if match_ind == args.hit_depth:
                break
            spectra = spectra_map[match["scan"]]
            const_mod_pos, const_mod_masses, n_variable = process_mods(
                args, match["mod_positions"], match["mod_masses"]
            )
            if n_variable > 0 and not np.isin(4294967295, const_mod_pos):
                if args.match_save:
                    save_match(spectra, match)
                ascore.score(
                    spectra["mz_values"], 
                    spectra["intensity_values"], 
                    match["peptide"], n_variable,
                    const_mod_pos, const_mod_masses
                )
                alt_sites = [",".join([str(site) for site in site_list]) 
                             for site_list in ascore.alt_sites]
                scores.append([match["scan"], 
                               ascore.best_sequence, 
                               ascore.best_score,
                               ";".join([str(s) for s in ascore.ascores]),
                               ";".join(alt_sites)])

    score_dataframe = DataFrame(scores,
                                columns=["Scan",
                                         "LocalizedSequence",
                                         "PepScore",
                                         "Ascores",
                                         "AltSites"])
    score_dataframe.to_csv(args.out_file, sep="\t", index=False)
    print("{} -- Ascore Completed".format(get_time_stamp()))


if __name__ == "__main__":
    main()

