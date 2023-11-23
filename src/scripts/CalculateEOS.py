import os
import numpy as np
import matplotlib.pyplot as plt


from ase.optimize import FIRE
from ase.db import connect
from ase.units import kJ
from ase.eos import EquationOfState
from ase.constraints import FixAtoms
from ase.constraints import StrainFilter, UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry


import copy


def to_scalar(value):
    if isinstance(value, np.ndarray) and value.size == 1:
        return value.item()
    else:
        return value


def CalculateEOS(structure,
             potential,
             of,
             dbname="NiTi_EOS.db",
             strain=np.linspace(-0.10,0.10,100),
             eos_eq='sj',
             fmax=2.0e-3,
             fstep=5000,
             ):
    """
    Calculate the equation of state for different potentials on a specified structure. 
    The structure is strained using the StrainFilter class from the ASE.
    
    TODO: Need to handle the sensitivities of the LAMMPSlib calc
    """
    db = connect(dbname)
#    db.serial = True
    potname = potential[0]
    structure.set_calculator(potential[1])
    vol_intial = structure.get_cell().volume

    # Store the symmetry
    sym_before = check_symmetry(structure,1.0e-6, verbose=False)
    # Apply SymmetryFilter and UnitCellFilter for optimization
    structure.set_constraint(FixSymmetry(structure,symprec=1.0e-4))
    ucf_sym = UnitCellFilter(structure)
    # FIRE is used as best performance seen.
    opt = FIRE(ucf_sym)
    opt.run(fmax=fmax,steps=fstep)
    #ISSUE: Need to turn of constraints for ASE db
    structure.set_constraint(None)
    etotal = structure.get_total_energy()
    stress = structure.get_stress()
    forces = structure.get_forces()
    ecoh = structure.get_potential_energy() / structure.get_global_number_of_atoms()
    nextrow = db.write(structure,
                       relaxed=True,
                       vi=vol_intial,
                       ecoh = ecoh,
                       spacegroup=sym_before["number"],
                       structure_name=structure.info['structure_name'],
                       model_name = structure.info['model_name'],
                       )

             

    sym_after = check_symmetry(structure,1.0e-6,verbose=False)

    if sym_before["number"] != sym_after["number"]:
        print("Not adding structure as the spacegroup symmetry has not been maintained!")
        return None
    
    # Now take optimized structure, remove symmetry constraint, copy, set calc, and calculate EOS
    structure_opt = structure.copy()
    #structure_opt.set_calculator(potential[1])
    sopt = structure_opt.copy()
    sopt.set_calculator(potential[1])
    sf = StrainFilter(sopt)
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
        
    # Plot and save EOS
    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_subplot(1, 1, 1)
    structure_name = structure.info["structure_name"]
    eos.plot(ax=ax,
             filename=os.path.join(of,f'{potname}_NiTi_{structure_name}_EOS_{eos_eq}.png'),
             )
    fig.clear()
    plt.close(fig)  

    
    db.update(id=nextrow,
              v0=v0,
              e0=e0,
              bulk_modulus=B / kJ * 1.0e24,
              fit_eq=eos_eq,
              data={"volumes":volumes,
                   "energies":energies},
             )
    
    # db.write(structure,
    #          relaxed=True,
    #          vi=vol_intial,
    #          v0=v0,
    #          e0=e0,
    #          ecoh=ecoh,
    #          etotal=etotal,
    #          bulk_modulus=B / kJ * 1.0e24,
    #          fit_eq=eos_eq,
    #          spacegroup=sym_before["number"],
    #          structure_name=structure.info['structure_name'],
    #          model_name = structure.info['model_name'],
    #          data={"volumes":volumes,
    #                "energies":energies},
    #          )

    return None


