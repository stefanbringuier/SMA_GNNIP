from ase.db import connect
from ase.optimize import FIRE
from ase.constraints import UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry

from ase.units import m, kg

grams = kg * (1 / 1000)
cm = m * (1 / 100)


def minimize_structure(structure, potential, write_db=None, fmax=5.0e-4, fstep=5000):
    """
    Minimize the ASE calculator using specific model and calculator.

    Arguments:
        structure (ASE.Atoms): structure to be minimized
        potential (Tuple(string,ASE.calculators)): Name of the potential and ASE calculator
        write_db (Union[None,string]): Optional, name of the ASE database to use.

    Returns:
        Union(int,(int,ASE.Atoms)): returns current row id or/and the ASE.Atoms structure

    """

    structure.set_calculator(potential[1])
    vol_intial = structure.get_cell().volume

    # Store the symmetry
    sym_before = check_symmetry(structure, 1.0e-6, verbose=False)

    # Apply SymmetryFilter and UnitCellFilter for optimization
    structure.set_constraint(FixSymmetry(structure, symprec=1.0e-4))
    ucf_sym = UnitCellFilter(structure)

    # FIRE is used as best performance seen.
    opt = FIRE(ucf_sym)
    opt.run(fmax=fmax, steps=fstep)

    sym_after = check_symmetry(structure, 1.0e-6, verbose=False)

    if sym_before["number"] != sym_after["number"]:
        ValueError(
            "Structure spacegroup symmetry has not been maintained after optimization!"
        )

    # NOTE: Need to turn of constraints for ASE db
    structure.set_constraint(None)

    if write_db:
        db = connect(write_db)
        chemsys = "".join(sorted(set(structure.symbols)))
        etotal = structure.get_total_energy()
        stress = structure.get_stress()
        forces = structure.get_forces()
        ecoh = structure.get_potential_energy() / structure.get_global_number_of_atoms()

        rho = structure.get_masses().sum() / structure.get_volume()
        rho *= cm**3 / grams  # gives g/cm^3

        nextrow = db.write(
            structure,
            relaxed=True,
            chemsys=chemsys,
            vol_unrelaxed=vol_intial,
            density=rho,
            ecoh=ecoh,
            spacegroup=sym_before["number"],
            structure_name=structure.info["structure_name"],
            model_name=structure.info["model_name"],
        )
        return nextrow, structure.copy()

    else:
        return structure.copy()
