import paths

from ase.db import connect
import ase.units as units
from ase.spacegroup import crystal


#from parcalc import ParCalculate
from elastic import get_elementary_deformations
from elastic.elastic import get_cij_order
from elastic import get_elastic_tensor

from ase.calculators.lammpslib import LAMMPSlib

from queue import Empty

from multiprocessing import Process, Queue

import time
import os
import tempfile
import shutil
from copy import deepcopy
from subprocess import check_output
from contextlib import contextmanager


def SerialCalculate(systems,calc,cleanup=True,prefix="Calc_"):

    if type(systems) != type([]) :
        sysl=[systems]
    else :
        sysl=systems
    basedir=os.getcwd()
    res=[]
    for n,s in enumerate(sysl):
        if isinstance(calc, LAMMPSlib):
            s.set_calculator(calc)
        else:
            s.set_calculator(deepcopy(calc))
        s.get_calculator().block=False
        place=tempfile.mkdtemp(prefix=prefix, dir=basedir)
        os.chdir(place)
        s.get_calculator().working_dir=place
        #print("Start at :", place)
        s.get_potential_energy()
        os.chdir(basedir)
        time.sleep(0.2)
        res.append([n,s])
        print("Workers started:", len(sysl))

    return [r for ns,s in enumerate(sysl) for nr,r in res if nr==ns]


def calculate_elastic_constants(
    dbname, structure, potential, npoints=10, displacement=2.0
):
    """
    displacement: percentage of length and angles
    npoints: number of iterations along displacement
    """
    
    db = connect(dbname)
    potname = potential[0]
    calc = potential[1]

        
        
    for entry in db.select(structure_name=structure, model_name=potname):
        spg = entry.spacegroup
        structure = crystal(entry.toatoms(), spacegroup=spg)
        # Need to set some defaults
        calc.clean = potential[1].reset
        # calc.set(isif=2)
        structure.calc = calc

        systems = get_elementary_deformations(structure, n=npoints, d=displacement)

        # Run the stress calculations on deformed cells
        calcpath = str(paths.data / potname.upper() / "Elastic_Calc_")
        result = SerialCalculate(systems, calc, prefix=calcpath)

        # Elastic tensor by internal routine
        ordering = get_cij_order(structure)
        cij, _ = get_elastic_tensor(structure, systems=result)
        cij /= units.GPa

        C = dict(zip(ordering, cij))
        db.update(entry.id, data={"elastic_constants": C})

    return None
