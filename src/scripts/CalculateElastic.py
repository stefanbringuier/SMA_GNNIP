import paths

from ase.db import connect
import ase.units as units
from ase.spacegroup import crystal


from parcalc import ParCalculate
from elastic import get_elementary_deformations
from elastic.elastic import get_cij_order
from elastic import get_elastic_tensor


def CalculateElasticConstants(dbname,
                  potential,
                  npoints=10,
                  displacement=0.33):

    db = connect(dbname)
    potname = potential[0]
    calc = potential[1]
    for entry in db.select(model_name = potname):
        spg = entry.spacegroup
        structure = crystal(entry.toatoms(),spacegroup=spg)
        # Need to set some defaults
        calc.clean = potential[1].reset
        #calc.set(isif=2)
        structure.calc = calc

        systems = get_elementary_deformations(structure, n=npoints, d=displacement)

        # Run the stress calculations on deformed cells
        calcpath = str(paths.data / potname.upper() / "Calc_")
        result = ParCalculate(systems,calc,block=False,prefix=calcpath)

        # Elastic tensor by internal routine
        ordering = get_cij_order(structure)
        cij, _ = get_elastic_tensor(structure, systems=result)
        cij /= units.GPa

        C = dict(zip(ordering,cij))
        db.update(entry.id,data={"elastic_constants" : C})

    return None
