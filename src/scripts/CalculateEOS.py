import os
import paths

import numpy as np
import matplotlib.pyplot as plt

from ase.db import connect
from ase.constraints import StrainFilter
from ase.units import kJ
from ase.eos import EquationOfState


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
    """Calculates the equation of state (EOS) for a given structure under strain.

    This function applies isotropic strain to the given structure, computes the energy
    and volume for each strain state, and fits these to an equation of state. Optionally,
    it plots the EOS curve.

    Args:
        structure (ASE.Atoms): The structure to be analyzed.
        potential (Tuple[str, ASE.calculators]): A tuple containing the potential type
            and its associated calculator.
        strain (List[float]): A list of strain values to be applied to the structure.
        eos_eq (str, optional): The type of EOS to be fitted (e.g., "sj" for Sjevold-Johnson).
            Defaults to "sj".
        plot (bool, optional): If True, the EOS curve will be plotted. Defaults to True.

    Returns:
        tuple:
            - v0 (float): Equilibrium volume.
            - e0 (float): Equilibrium energy.
            - B (float): Bulk modulus.
            - volumes (List[float]): List of volumes corresponding to the strains.
            - energies (List[float]): List of energies corresponding to the strains.
    """

    
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
    """Minimizes the given structure, calculates its equation of state (EOS), and updates a database.

    This function first minimizes the provided structure using the specified potential, then calculates
    the EOS based on the optimized structure and a given strain range. The resulting EOS parameters and
    data are updated in the specified database.

    Args:
        structure (ASE.Atoms): The initial structure to be minimized and analyzed.
        potential (Tuple[str,ASE.calculators]): The potential used for minimizing the structure.
        dbname (str, optional): The name of the database file to store results. Defaults to "EOS.db".
        strain (np.ndarray, optional): An array of strain values for which the EOS will be calculated.
            Defaults to a linear space between -0.10 and 0.10 with 40 points.

    Returns:
        None: This function does not return a value but updates the database with EOS data.
    """



    nextrow, opt_structure = minimize_structure(structure,potential,write_db=dbname)
             
    v0,e0,B,volumes,energies = calculate_eos(opt_structure,potential,strain)

    db = connect(dbname)
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
