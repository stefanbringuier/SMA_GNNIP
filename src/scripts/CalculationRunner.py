import warnings
import argparse
import os
import paths

import numpy as np

from ase.db import connect
from Structures import get_structure
from Calculators import get_ase_calculator
from CalculateEOS import get_min_structure_and_calculate_eos
from CalculateElastic import calculate_elastic_constants
from CalculatePhonons import calculate_phonons

# Surpress warning from dependencies for which I can do nothing about.
warnings.filterwarnings("ignore", category=UserWarning, module="ase.spacegroup")
warnings.filterwarnings(
    "ignore", category=UserWarning, module="dgl.backend.pytorch.tensor"
)
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_cuda'.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_cudnn'.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_mps'.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_mkldnn'.*")


def min_eos_calc(dbname, structure, potential, strain_info):
    # NOTE: Because this adds a new entry and its hard to get Snakemake to handle
    # checking if the entry exist directly. We check here and then skip if the model
    # and structure have already been optimized
    entry_exist = False
    try:
        entry = connect(dbname).get(model_name=potential[0], structure_name=structure)
        if entry != None:
            entry_exist = True
    except:
        print(
            f"""
            -----------------------------------------------------------
            Started EOS Calculation of {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
            -----------------------------------------------------------
            """
        )
        get_min_structure_and_calculate_eos(
            structure, potential, dbname=dbname, strain=np.linspace(*strain_info)
        )
        print(
            f"""
            ------------------------------------------------------------
            Finished EOS Calculation of {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
            ------------------------------------------------------------
            """
        )
    return entry_exist


def phonon_calc(dbname, structure, potential, strain_info):
    """
    dbname (string_info)
    strain Tuple(min,max,number)
    potential Tuple(model_name, asecalc)
    """
    print(
        f"""
    -----------------------------------------------------------------
    Started Phonon Calculation of {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
    -----------------------------------------------------------------
    """
    )
    strain = np.linspace(*strain_info)
    calculate_phonons(
        dbname, structure.info["structure_name"], potential, strain=strain
    )
    print(
        f"""
    ------------------------------------------------------------------
    Finished Phonon Calculation of {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
    ------------------------------------------------------------------
    """
    )

    return None


def elastic_calc(dbname, structure, potential):
    print(
        f"""
    ----------------------------------------------------------------------------
    Started Elastic Constants Calculation of {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
    ----------------------------------------------------------------------------
    """
    )
    # NOT FULLY VALIDATED
    calculate_elastic_constants(dbname, structure.info["structure_name"], potential)
    print(
        f"""
    -----------------------------------------------------------------------------
    Finished Elastic Constants Calculation of {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
    -----------------------------------------------------------------------------
    """
    )

    return None


def simulation():
    return ValueError("Simulation routine not implemented!")


def main():
    """
    Provides equation of state calculations of material. Be mindful that correct python environment is
    activated for selected model (i.e., ASE Calculator) otherwise routines will fail to run.
    """
    parser = argparse.ArgumentParser(
        description="""Top-level run script for EOS, Phonons, and Elastic constants.
    Be mindful that the correct python environment is activated for the model (i.e., ASE Calculator) you
    intend to use!!"""
    )
    parser.add_argument(
        "--dbname", type=str, default="Results.json", help="ASE Database file name"
    )
    parser.add_argument("--dbfolder", type=str, default=".", help="ASE Database folder")
    parser.add_argument(
        "--model",
        choices=[
            "Mutter",
            "MutterASE",
            "Zhong",
            "ZhongASE",
            "Ko",
            "M3GNet",
            "CHGNet",
            "MACE",
            "ALIGNN",
        ],
        default="Mutter",
        help="Choose ASECalculator Model",
    )
    parser.add_argument("--structure", default="B2")
    parser.add_argument(
        "--n",
        type=int,
        default=40,
        help="Number of configurations in strain scan of EOS",
    )
    parser.add_argument(
        "--nph",
        type=int,
        default=13,
        help="Number of configurations in strain scan of Phonons",
    )
    parser.add_argument(
        "--min_strain_eos",
        type=float,
        default=-0.10,
        help="Minimum strain for eos calcs, fractional",
    )
    parser.add_argument(
        "--max_strain_eos",
        type=float,
        default=0.10,
        help="Maximum strain for eos calcs, fractional",
    )
    parser.add_argument(
        "--min_strain_ph",
        type=float,
        default=-0.02,
        help="Minimum strain for phonon calcs, fractional",
    )
    parser.add_argument(
        "--max_strain_ph",
        type=float,
        default=0.02,
        help="Maximum strain for phonon calcs, fractional",
    )
    parser.add_argument("--eos_fit", type=str, default="sj", help="EOS fit function")
    parser.add_argument(
        "--calc_type",
        type=str,
        choices=["min_eos", "phonons", "elastic", "simulation"],
        required=True,
        help="Specify the type of calculation to perform (required).",
    )

    args = parser.parse_args()

    outfolder = paths.data / args.model.upper()
    os.makedirs(outfolder, exist_ok=True)

    structure = get_structure(args.structure)
    asecalc = get_ase_calculator(args.model)
    structure.calc = asecalc
    structure.info["model_name"] = args.model

    db_path_file = paths.data / args.dbname

    if args.calc_type == "min_eos":
        min_eos_calc(
            db_path_file,
            structure,
            (args.model, asecalc),
            (args.min_strain_eos, args.max_strain_eos, args.n),
        )
    elif args.calc_type == "phonons":
        phonon_calc(
            db_path_file,
            structure,
            (args.model, asecalc),
            (args.min_strain_ph, args.max_strain_ph, args.nph),
        )
    elif args.calc_type == "elastic":
        elastic_calc(db_path_file, structure, (args.model, asecalc))

    return None


if __name__ == "__main__":
    main()
