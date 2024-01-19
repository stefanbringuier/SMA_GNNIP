import csv
from copy import deepcopy

import ase.units as units
import paths
from ase.calculators.lammpslib import LAMMPSlib
from ase.db import connect
from ase.spacegroup import crystal
from deepmd.calculator import DP
from elastic import get_elastic_tensor, get_elementary_deformations
from elastic.elastic import get_cij_order


def serial_calculate(systems, calculator):
    """
    Create the list of systems with a calculator for elastic constant calculations. This is done in serial.

    Args:
        systems (Union[list(ASE.Atoms),ASE.Atoms]): atomic structures that contain deformations
        calculator (ASE.calculators): Internal ASE calculator instance

    Returns:
        list(systems)
    """

    if type(systems) != type([]):
        sysl = [systems]
    else:
        sysl = systems

    res = []
    for n, s in enumerate(sysl):
        if isinstance(calculator, LAMMPSlib):
            s.set_calculator(calculator)
        elif isinstance(calculator, DP):
            s.set_calculator(calculator)
        else:
            s.set_calculator(deepcopy(calculator))

        s.get_potential_energy()
        res.append([n, s])

    return [r for ns, s in enumerate(sysl) for nr, r in res if nr == ns]


def calculate_elastic_constants(
    dbname,
    chemsys,
    structure_name,
    potential,
    npoints=10,
    displacement=0.25,
    update=True,
    write_file=True,
):
    """
    Calculate the elastic constants for a given structure and potential and write it to an
    ASE database file. The structure is retrieved from the database.

    Args:
        dbname (string): The path of the ASE database file.
        chemsys (string): Chemical system, e.g., NiTi, PtTi, NiAl3
        structure_name (string): The name of the structure in the ASE database file.
        potential (string): The name of the potential style.
        npoints (int,10):  Optional, sets the number of displacement points.
        displacement (int/float,0.5): Optional, the min and max atomic displacement in percent.
        update (bool,True): Optional, wheter or not to update ASE database file.
        write_file (bool,True): Write elastic constants to csv

    Returns:
        None

    Notes:
        Writes the results from the calculation to `dbname`
    """

    db = connect(dbname)
    potname = potential[0]
    calculator = potential[1]

    entry = db.get(
        chemsys=chemsys, structure_name=structure_name, model_name=potname, relaxed=True
    )

    spg = entry.spacegroup
    structure = crystal(entry.toatoms(), spacegroup=spg)
    # calculator.clean = potential[1].reset
    structure.calc = calculator

    systems = get_elementary_deformations(structure, n=npoints, d=displacement)

    result = serial_calculate(systems, calculator)

    ordering = get_cij_order(structure)
    cij, _ = get_elastic_tensor(structure, systems=result)
    cij /= units.GPa

    C = dict(zip(ordering, cij))

    if write_file:
        filename = str(
            paths.data
            / potname.upper()
            / f"{chemsys}_{structure_name}_{potname}_Cij.csv"
        )
        with open(filename, mode="w", newline="") as file:
            writer = csv.writer(file)
            writer.writerow(["Elastic Constant", "Value [GPa]"])
            for key, value in C.items():
                writer.writerow([key, value])
    if update:
        db.update(entry.id, data={"elastic_constants": C})
    else:
        return C

    return None


if __name__ == "__main__":
    # Testing data
    print("!!! TESTING ElasticCalculation.py !!!")
    from Calculators import *

    calculator = get_ase_calculator("M3GNet")
    dbname = paths.data / "Results.json"
    structure_name = "BCO"
    chemsys = "NiTi"
    potential = ("M3GNet", calculator)
    for d in [0.25, 0.5]:
        for n in [10, 15]:
            C = calculate_elastic_constants(
                dbname,
                chemsys,
                structure_name,
                potential,
                npoints=n,
                displacement=d,
                update=False,
            )
            print(f"{d} {n} {C}")
    #    assert list(C.values()) == [0.0,0.0,0.0]
