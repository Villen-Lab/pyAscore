import sys
import argparse
import re
from datetime import datetime
from itertools import groupby
import numpy as np
from pandas import DataFrame
from tqdm import tqdm
import pickle
from pyascore import *
from .config import *

def get_time_stamp():
    return datetime.now().strftime("%m/%d/%y %H:%M:%S")

def parse_spectra(arg_ref):
    print("{} -- Reading spectra from: {}".format(get_time_stamp(), arg_ref.spec_file))
    spec_parser = SpectraParser(arg_ref.spec_file, arg_ref.spec_file_type)
    return spec_parser.to_dict()


def parse_identifications(arg_ref):
    print("{} -- Reading identifications from: {}".format(get_time_stamp(), arg_ref.ident_file))

    # Update static modifications
    static_mods = {}
    for aa_group, mass in zip(arg_ref.static_mod_groups.split(","),
                              arg_ref.static_mod_masses.split(",")):
        static_mods.update({aa : float(mass) for aa in aa_group})

    # Build mass corrector
    mods = COMMON_MODS.copy()
    mods.update({aa : arg_ref.mod_mass for aa in arg_ref.residues})
    mods.update(static_mods)
    mass_corrector = MassCorrector(mod_mass_dict=mods)

    # Parse
    id_parser = iter(
        sorted(IdentificationParser(arg_ref.ident_file,
                                    arg_ref.ident_file_type,
                                    mass_corrector,
                                    static_mods=static_mods).to_list(),
               key=lambda spec: spec["scan"])
    )
    return id_parser


def validate_args(arg_ref):
    # Check modifiable resiudes
    allowed_residues = "ncACDEFGHIKLMNOPQRSTUVWY"
    for aa in arg_ref.residues:
        if aa not in allowed_residues:
            raise ValueError("The residue inputed, {}, is not allowed."
                             " Must be one of: {}".format(aa, allowed_residues))

    # Check fragment types
    allowed_fragments = "cbyzZ"
    for frag in arg_ref.fragment_types:
        if frag not in allowed_fragments:
            raise ValueError("The fragment type inputed, {}, is not allowed."
                             " Must be one of: {}".format(frag, allowed_fragments))

    # Check max fragment charge
    if arg_ref.max_fragment_charge < 1:
        raise ValueError("The max fragment charge must be greater than or equal to 1")


def build_ascore(arg_ref):
    ascore = PyAscore(bin_size=100., n_top=10,
                      mod_group=arg_ref.residues,
                      mod_mass=arg_ref.mod_mass,
                      mz_error=arg_ref.mz_error,
                      fragment_types=arg_ref.fragment_types)

    if arg_ref.neutral_loss_masses and arg_ref.neutral_loss_masses:
        nl_groups = arg_ref.neutral_loss_groups.split(",")
        nl_masses = [float(m) for m in arg_ref.neutral_loss_masses.split(",")]
        [ascore.add_neutral_loss(g, m) for g,m in zip(nl_groups, nl_masses)]
 
    return ascore


def process_mods(arg_ref, sequence, positions, masses):
    variable_mod_count = 0
    const_pos, const_masses = [], []
    for pos, mass in zip(positions, masses):
        shift = 1 if arg_ref.zero_based else 0
        if pos + shift == 0:
            aa = "n"
        else:
            aa = sequence[pos - 1 + shift]

        matches_variable = np.isclose(arg_ref.mod_mass, mass,
                                      rtol=1e-6, atol=arg_ref.mod_correction_tol)
        if matches_variable and aa in arg_ref.residues:
            variable_mod_count += 1
        else:
            const_pos.append(pos + shift)
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
    parser = build_parser()
    args = parser.parse_args()
    if args.parameter_file:
        args = parser.parse_args(args_from_file(args.parameter_file) + sys.argv[1:])
    validate_args(args)

    print("{} -- Ascore Started".format(get_time_stamp()))

    spectra_map = parse_spectra(args)
    psms = parse_identifications(args)

    print("{} -- Anlyzing PSMs".format(get_time_stamp()))
    ascore = build_ascore(args)
    scores = []
    for _, match_group in tqdm(groupby(psms, lambda x: x["scan"])):
        for match_ind, match in enumerate(match_group):
            if match_ind == args.hit_depth:
                break
            spectra = spectra_map[match["scan"]]
            const_mod_pos, const_mod_masses, n_variable = process_mods(
                args, match["peptide"], match["mod_positions"], match["mod_masses"]
            )

            # Try and figure out PSM charge
            if match["charge_state"] is not None and match["charge_state"] != 0:
                psm_charge = match["charge_state"]
            elif spectra["precursor_charge"] is not None and spectra["precursor_charge"] != 0:
                psm_charge = spectra["precursor_charge"]
            else:
                psm_charge = 2
            psm_charge = max(psm_charge, 2)

            if n_variable > 0:
                if args.match_save:
                    save_match(spectra, match)
                ascore.score(
                    spectra["mz_values"], 
                    spectra["intensity_values"], 
                    match["peptide"], n_variable,
                    min(args.max_fragment_charge,
                        psm_charge - 1),
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

