import paths
from copy import deepcopy

from ase.db import connect
import ase.units as units
from ase.spacegroup import crystal
from ase.calculators.lammpslib import LAMMPSlib

from elastic import get_elementary_deformations
from elastic.elastic import get_cij_order
from elastic import get_elastic_tensor


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
        else:
            s.set_calculator(deepcopy(calculator))

        s.get_potential_energy()
        res.append([n, s])

    return [r for ns, s in enumerate(sysl) for nr, r in res if nr == ns]


def calculate_elastic_constants(
    dbname,
    structure_name,
    potential,
    npoints=10,
    displacement=0.5,
    update=True,
):
    """
    Calculate the elastic constants for a given structure and potential and write it to an
    ASE database file. The structure is retrieved from the database.

    Args:
        dbname (string): The path of the ASE database file.
        structure_name (string): The name of the structure in the ASE database file.
        potential (string): The name of the potential style.
        npoints (int,10):  Optional, sets the number of displacement points.
        displacement (int/float,0.5): Optional, the min and max atomic displacement in percent.
        update (bool,True): Optional, wheter or not to update ASE database file.

    Returns:
        None

    Notes:
        Writes the results from the calculation to `dbname`
    """

    db = connect(dbname)
    potname = potential[0]
    calculator = potential[1]

    for entry in db.select(structure_name=structure_name, model_name=potname):
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
    dbname = paths.data / "NiTi_Structures.json"
    structure_name = "BCO"
    potential = ("M3GNet",calculator)
    for d in [0.25,0.5]:
        for n in [10,15]:
            C = calculate_elastic_constants(dbname, structure_name,potential,npoints=n,displacement=d,update=False)
            print(f"{d} {n} {C}")
    #    assert list(C.values()) == [0.0,0.0,0.0]
