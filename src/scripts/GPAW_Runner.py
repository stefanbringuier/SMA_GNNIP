import argparse
import os
import sys

sys.path.append(os.path.dirname(os.path.abspath(__file__)))

# Surpress warning from dependencies for which I can do nothing about.
import warnings

import numpy as np
import paths
from ase import Atoms
from ase.db import connect
from CalculateElastic import calculate_elastic_constants
from CalculateEOS import get_min_structure_and_calculate_eos
from CalculatePhonons import calculate_phonons
from Calculators import get_gpaw_calculator
from gpaw import mpi
from Structures import get_structure

warnings.filterwarnings("ignore", category=UserWarning, module="ase.spacegroup")


def isolated_atom(element, a=10.0):
    atom = Atoms(
        element,
        scaled_positions=[[0.5, 0.5, 0.5]],
        cell=[a, a, a],
        pbc=[True, True, True],
    )

    return atom


def min_eos_calc(dbname, structure, potential, strain_info):
    # NOTE: Because this adds a new entry and its hard to get Snakemake to handle
    # checking if the entry exist directly. We check here and then skip if the model
    # and structure have already been optimized
    chemsys = structure.info["chemsys"]
    structure_name = structure.info["structure_name"]
    entry_existed = False
    try:
        entry = connect(dbname).get(
            chemsys=chemsys,
            model_name=potential[0],
            structure_name=structure_name,
            relaxed=True,
        )
        print(
            f"""
            -----------------------------------------------------------
            EOS already exist for {structure.info['chemsys']} {structure.info['structure_name']} using {potential[0]}
            -----------------------------------------------------------
            """
        )
        entry_existed = True

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

    return entry_existed


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

    calculate_phonons(
        dbname,
        structure.info["chemsys"],
        structure.info["structure_name"],
        potential,
        strain=strain_info,
        plotting=True,
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
    calculate_elastic_constants(
        dbname, structure.info["chemsys"], structure.info["structure_name"], potential
    )
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="""Top-level run script for EOS, Phonons, and Elastic constants.
    Be mindful that the correct python environment is activated for GPAW
    """
    )
    parser.add_argument(
        "--dbname", type=str, default="Results.json", help="ASE Database file name"
    )
    parser.add_argument("--dbfolder", type=str, default=".", help="ASE Database folder")
    parser.add_argument("--structure", default="B2")
    parser.add_argument(
        "--n",
        type=int,
        default=20,
        help="Number of configurations in strain scan of EOS",
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
    parser.add_argument("--eos_fit", type=str, default="sj", help="EOS fit function")
    parser.add_argument("--chemsys", type=str, default="NiTi", help="Chemical system")

    args = parser.parse_args()
    outfolder = paths.data / "GPAW"
    os.makedirs(outfolder, exist_ok=True)

    structure = get_structure(args.chemsys, args.structure)

    # # Get isolated atom energie
    # atoms = set(structure.get_chemical_symbols())
    # atom_energies = {}
    # gpawcalc = get_gpaw_calculator()#kpts=(1,1,1))
    # for a in atoms:
    #     atom = isolated_atom(a)
    #     atom.set_calculator(gpawcalc)
    #     atom_energies[a] = atom.get_total_energy()
    # structure.info["dft_atom_energies"] = atom_energies

    gpawcalc = get_gpaw_calculator()
    structure.calc = gpawcalc
    structure.info["model_name"] = "GPAW"

    db_path_file = paths.data / args.dbname

    min_eos_calc(
        db_path_file,
        structure,
        ("GPAW", gpawcalc),
        (args.min_strain_eos, args.max_strain_eos, args.n),
    )

    # # Connect to db and shift energies based on atom energies
    # db = connect(db_path_file)
    # entry = db.get(
    #     model_name="GPAW", structure_name=args.structure, chemsys=args.chemsys
    # )
    # ecoh = entry.ecoh - sum(atom_energies.values())
    # eos_energies = [e - sum(atom_energies.values()) for e in entry.data["energies"]]
    # db.update(
    #     id=entry.id,
    #     ecoh=ecoh,
    #     data={"energies": eos_energies,
    #           "atom_energies":atom_energies,
    #           }
    # )

    # We only run 0.0 strain case
    gpawcalc = get_gpaw_calculator(kpts=(4, 4, 4))
    phonon_calc(
        db_path_file,
        structure,
        ("GPAW", gpawcalc),
        [0.00],
    )

    elastic_calc(db_path_file, structure, ("GPAW", gpawcalc))
