import os
import paths

import numpy as np
import matplotlib.pyplot as plt

from ase.db import connect
from ase.constraints import StrainFilter
from ase.units import kJ
from ase.eos import EquationOfState

import copy


def to_scalar(value):
    if isinstance(value, np.ndarray) and value.size == 1:
        return value.item()
    else:
        return value

def calculate_eos(structure,
                  potential,
                  strain,
                  eos_eq="sj",
                  plot=True,
                  ):
    
    structure.set_calculator(potential[1])

    sf = StrainFilter(structure)
    energies = []
    volumes = []
    for e in strain:
        isotropic_strain = np.array([e,e,e, 0, 0, 0])
        sf.set_positions(isotropic_strain)
        energy = to_scalar(sf.get_potential_energy())
        volume = sf.get_cell().volume
        energies.append(energy)
        volumes.append(volume)
        
    eos = EquationOfState(volumes, energies,eos=eos_eq)
    v0, e0, B = eos.fit()


    if plot:
        species = "".join(set(structure.get_chemical_symbols()))
        name = structure.info["structure_name"]
        model = potential[0]
        plot_eos(eos,model,species,name)
        

    return (v0, e0, B, volumes, energies)

def plot_eos(eos,potential_name,specie_names,structure_name):
    # Plot and save EOS
    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_subplot(1, 1, 1)
    file_name = f"{specie_names}_{structure_name}_{potential_name}_EOS.png"
    file_path = os.path.join(potential_name.upper(),file_name)
    plot_file = str(paths.data / file_path)
    eos.plot(ax=ax,filename=plot_file)
    fig.clear()
    plt.close(fig)  
    
def get_min_structure_and_calculate_eos(structure,
             potential,
             dbname="EOS.db",
             strain=np.linspace(-0.10,0.10,40),
             ):
    """

    """

    nextrow, opt_structure = minimize_structure(structure,potential,write_db=dbname)
             
    v0,e0,B,volumes,energies = calculate_eos(opt_structure,potential,strain)
    
    db.update(id=nextrow,
              v0=v0,
              e0=e0,
              bulk_modulus=B / kJ * 1.0e24,
              fit_eq=eos_eq,
              data={"volumes":volumes,
                   "energies":energies},
             )

    
    return None


if __name__ == "__main__":
    from Structures import NiTi_B2_Structure
    from Calculators import get_ase_calculator
    from MinimizeStructure import minimize_structure
    B2 = NiTi_B2_Structure()
    model = "MutterASE"
    calculator = get_ase_calculator(model)
    potential = (model,calculator)
    opt_structure = minimize_structure(B2,potential)
    strain = np.linspace(-0.10,0.10,100)
    v0, e0, B, v, e = calculate_eos(opt_structure,potential,strain)
