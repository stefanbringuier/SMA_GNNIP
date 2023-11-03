import argparse
import os
import numpy as np
import matplotlib.pyplot as plt
#from pathlib import Path
import paths



from ase.optimize import BFGS, FIRE, GPMin
from ase.db import connect
from ase.spacegroup import crystal
from ase.calculators.eam import EAM
from ase.units import kJ
from ase.eos import EquationOfState
from ase.constraints import StrainFilter, UnitCellFilter
from ase.spacegroup.symmetrize import FixSymmetry, check_symmetry

from NiTiStructures import *
from Calculators import *


def calc_EOS(structure,
             potential,
             of,
             dbname="NiTi_EOS.db",
             strain=np.linspace(-0.10,0.10,100),
             eos_eq='sj',
             fmax=5.0e-3,
             ):
    """
    Calculate the equation of state for different potentials on a specified structure. 
    The structure is strained using the StrainFilter class from the ASE.
    """


    potname = potential[0]
    structure.calc =  potential[1]
    vol_intial = structure.get_cell().volume
    
    db = connect(dbname)

    sf = StrainFilter(structure)
    energies = []
    volumes = []
    for e in strain:
        isotropic_strain = np.array([e,e,e, 0, 0, 0])
        sf.set_positions(isotropic_strain)
        energy = sf.get_potential_energy()
        volume = sf.get_cell().volume
        energies.append(energy)
        volumes.append(volume)
        
    eos = EquationOfState(volumes, energies,eos=eos_eq)
    v0, e0, B = eos.fit()


        
    # Plot and save EOS
    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_subplot(1, 1, 1)
    structure_name = structure.info["structure_name"]
    eos.plot(ax=ax, filename=os.path.join(of,f'{potname}_NiTi_{structure_name}_EOS_{eos_eq}.png'))
    fig.clear()
    plt.close(fig)  



    # Now find the ground-state structure using the StrainFilter and BFGS optimization
    # Then store the ground-state optimal structure and the results for the EOS fit.
    structure_copy = structure.copy()
    structure_copy.calc = potential[1]

    # Apply SymmetryFilter and UnitCellFilter for optimization
    structure_copy.set_constraint(FixSymmetry(structure_copy,symprec=1.0e-3))
    ucf_sym = UnitCellFilter(structure_copy)
    # Set up the BFGS optimizer
    opt = FIRE(ucf_sym)
    opt.run(fmax=fmax,steps=10000)
    ecoh = structure_copy.get_potential_energy() / structure_copy.get_global_number_of_atoms()

    # We print out the initial symmetry groups at two different precision levels
    #optimized_structure = ucf_sym.get_atoms()
    sym_before = check_symmetry(structure,1.0e-6, verbose=False)
    sym_after = check_symmetry(structure_copy,1.0e-6,verbose=False)
    if sym_before["number"] != sym_after["number"]:
        print("Not adding structure as the spacegroup symmetry has not been maintained!")
    else:
        # ISSUE: Need to turn of constraints for db
        structure_copy.set_constraint(None)
        db.write(structure_copy,
                 relaxed=True,
                 vi=vol_intial,
                 v0=v0,
                 e0=e0,
                 ecoh=ecoh,
                 bulk_modulus=B / kJ * 1.0e24,
                 fit_eq=eos_eq,
                 spacegroup=sym_before["number"],
                 structure_name=structure.info['structure_name'],
                 model_name = structure.info['model_name'],
                 data={"volumes":volumes,
                       "energies":energies},
                 )

    return None


def main():
    parser = argparse.ArgumentParser(description="""Equation of State Calculations of NiTi.
    Be mindful that the correct python environment is activated for the model (i.e., ASE Calculator) you
    intend to use!!""")
    parser.add_argument('--dbname', type=str, default='NiTi_EOS.json', help='ASE Database file name')
    parser.add_argument('--dbfolder', type=str, default='.', help='ASE Database folder')
    parser.add_argument('--model', choices=['Mutter', 'Zhong', 'M3GNet', 'CHGNet','MACE','ALIGNN'], default='Mutter', help='Choose ASECalculator Model')
    parser.add_argument('--structure',default='B2')
    parser.add_argument('--n',type=int,default=25,help='Number of configurations in strain scan')
    parser.add_argument('--min_strain', type=float, default=-0.1, help='Minimum strain, fractional')
    parser.add_argument('--max_strain', type=float, default=0.1, help='Maximum strain, fractional')
    parser.add_argument('--eos_fit', type=str, default='sj', help='EOS fit function')
    args = parser.parse_args()
    # Inside src/scripts/GenerateEOS.py



    banner = '''\
    ********************************************************
    * Equation of State Calculations of NiTi              *
    * Be mindful that correct python environment is       *
    * activated for selected model (i.e., ASE Calculator) *
    * For help: GenerateEOS.py -h                         *
    ********************************************************
    '''
    print(banner)


    outfolder = paths.data / args.model.upper()
    os.makedirs(outfolder, exist_ok=True)

    asecalc = GetCalculator(args)
    structure = GetStructure(args.structure)
    structure.calc = asecalc
    structure.info["model_name"] = args.model
    
    
    db_path_file = paths.data / args.dbname
    calc_EOS(structure,
             (args.model, asecalc),
             outfolder,
             strain=np.linspace(args.min_strain,args.max_strain,args.n),
             dbname=db_path_file,
             eos_eq=args.eos_fit)


if __name__ == "__main__":
    main()
