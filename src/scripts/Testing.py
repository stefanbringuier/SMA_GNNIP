import warnings
import argparse
import os
import paths

from NiTiStructures import *
from Calculators import *

from CalculateEOS import *

#from CalculateElastic import *
from CalculatePhonons import *

# Surpress warning from dependencies for which I can do nothing about.
warnings.filterwarnings("ignore", category=UserWarning, module="ase.spacegroup")
warnings.filterwarnings(
    "ignore", category=UserWarning, module="dgl.backend.pytorch.tensor"
)
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_cuda'.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_cudnn'.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_mps'.*")
warnings.filterwarnings("ignore", category=UserWarning, message=".*'has_mkldnn'.*")


def main():
    """
    Provides equation of state calculations of NiTi. Be mindful that correct python environment is
    activated for selected model (i.e., ASE Calculator)
    """
    parser = argparse.ArgumentParser(
        description="""Equation of State Calculations of NiTi.
    Be mindful that the correct python environment is activated for the model (i.e., ASE Calculator) you
    intend to use!!"""
    )
    parser.add_argument(
        "--dbname", type=str, default="NiTi_EOS.json", help="ASE Database file name"
    )
    parser.add_argument("--dbfolder", type=str, default=".", help="ASE Database folder")
    parser.add_argument(
        "--model",
        choices=["Mutter", "Zhong", "M3GNet", "CHGNet", "MACE", "ALIGNN"],
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
    args = parser.parse_args()

    outfolder = paths.data / args.model.upper()
    os.makedirs(outfolder, exist_ok=True)

    asecalc = GetCalculator(args)
    structure = GetStructure(args.structure)
    structure.calc = asecalc
    structure.info["model_name"] = args.model

    db_path_file = paths.data / args.dbname


    # NOTE: Because this adds a new entry and its hard to get Snakemake to handle
    # checking if the entry exist directly. We check here and then skip if the model
    # and structure have already been optimized
    entry_exist = False
    from ase.db import connect
    try:
        entry = connect(db_path_file).get(model_name=args.model, structure_name=args.structure)
        if entry != None:
            entry_exist = True 
    except:
        print(
            f"""
            -----------------------------------------------------
            EOS Calculation of {args.structure} using {args.model}
            -----------------------------------------------------
            """
        )
        CalculateEOS(
            structure,
            (args.model, asecalc),
            outfolder,
            strain=np.linspace(args.min_strain_eos, args.max_strain_eos, args.n),
            dbname=db_path_file,
            eos_eq=args.eos_fit,
        )

    # NOTE: Due to resource demands only B2 Phonons of ALIGNN
    # if args.model == "ALIGNN":
    #    if args.structure != "B2":
    #        return None

    print(
        f"""
    -----------------------------------------------------
    Phonon Calculation of {args.structure} using {args.model}
    -----------------------------------------------------
    """
    )
    strain = np.linspace(args.min_strain_ph, args.max_strain_ph, args.nph)
    CalculatePhonons(db_path_file, args.structure, (args.model, asecalc), strain=strain)

    print(
        f"""
    --------------------------------------------------------------------
    Elastic Constants Calculation of {args.structure} using {args.model}
    --------------------------------------------------------------------
    """
    )
    # NOT VALIDATED
#    calculate_elastic_constants(db_path_file, args.structure, (args.model, asecalc))

    return None


if __name__ == "__main__":
    main()
