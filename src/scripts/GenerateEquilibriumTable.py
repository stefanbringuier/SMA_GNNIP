import sys
import re
import numpy as np
from math import gcd
from functools import reduce
import paths


from ase.db import connect
from ase.spacegroup import get_basis
from ase.spacegroup import Spacegroup, get_spacegroup

from Config import *

def get_element_basis_list(structure, spacegroup, tol=1e-3):
    """Provide the basis for a crystal structure.

    Attempts to find the basis atoms for a given crystal structure. 

    STILL IN TESTING.

    Arguments:
       - structure (ASE.Atoms): The crystal structure.
       - spacegroup (int): The spacegroup number
       - tol (float,optional): The tolerance used for determining if sites are the same.

    returns:
         np.array: The basis atoms
    
    """
    sg = Spacegroup(spacegroup)
    rotations = sg.rotations
    translations = sg.translations

    unique_sites = []
    unique_elements = []

    for atom in structure:
        pos = atom.scaled_position
        is_unique = True

        # Apply symmetry operations
        for rot, trans in zip(rotations, translations):
            sym_pos = np.dot(rot, pos) + trans
            sym_pos -= np.floor(sym_pos)  # Wrap around using periodic boundary conditions

            # Check if the position is already in the unique sites within a tolerance
            for usite in unique_sites:
                if np.allclose(sym_pos, usite, atol=tol):
                    is_unique = False
                    break

            if not is_unique:
                break

        if is_unique:
            unique_sites.append(pos)
            unique_elements.append(atom.symbol)

    return list(zip(unique_elements, unique_sites))

def generate_structure_table(dbname, chemsys, output_filename):
    db = connect(dbname)
    
    unique_structure_types = set(entry.structure_name for entry in db.select())
    
    with open(output_filename, 'w') as file:
        file.write("\\begin{longtable}{|l|c|}\n")
        file.write("\\caption{Equilibrium structures for NiTi.} \\label{tab:equil_niti}\\\\\n")
        file.write("\\hline\n")
        file.write("Structure & Unit Cell \\\\\n")
        file.write("\\hline\n")
        file.write("\\endfirsthead\n")
        file.write("\\caption[]{Equilibrium structures for NiTi. (Continued)}\\\\\n")
        file.write("\\hline\n")
        file.write("Structure\n(Spacegroup \#) & Unit Cell \\\\\n")
        file.write("\\endhead\n")
        file.write("\\hline\n")
        file.write("\\endfoot\n")

        for structure_type in unique_structure_types:
            structure_name = structure_type.replace("_","-")
            spg_num = SPACEGROUP_MAP[structure_name]
            file.write(f"{structure_name}\n({spg_num}) & ")
            file.write("\\begin{tabular}{c c c c c c c}\n")  # Inner table header
            file.write("Model & $a$ (\AA) & $b$ (\AA) & $c$ (\AA) & $\\alpha$ ($^\circ$) & $\\beta$ ($^\circ$) & $\\gamma$ ($^\circ$)\\\\\n")

            entries = db.select(chemsys=chemsys, structure_name=structure_type)
            sentries = sorted(entries, key=lambda entry: ORDER[entry.model_name])
            for entry in sentries:
                structure = entry.toatoms()
                # Extract unit cell parameters
                a, b, c, alpha, beta, gamma = structure.cell.cellpar()
                model_name =  entry.get('model_name', 'NA')
                file.write("\\hline\\noalign{\\smallskip}\n")
                file.write(f"{model_name} & {a:.3f} & {b:.3f} & {c:.3f} & {alpha:.2f} & {beta:.2f} & {gamma:.2f}\\\\\n")
                # Atom positions
                file.write("& & & \\multicolumn{4}{c}{\\begin{tabular}{c c c c}\n")  # Corrected inner table header
                file.write("& x & y & z\\\\\n")
                file.write("\\hline\n")
                basis = get_element_basis_list(structure,entry.spacegroup)
                for e, b in basis:
                    file.write(f"{e} & {b[0]:.4f} & {b[1]:.4f} & {b[2]:.4f} \\\\\n")
                file.write("\\end{tabular}}\\\\\n")

            file.write("\\end{tabular} \\\\\n")
            file.write("\\hline\n")

        file.write("\\end{longtable}\n")
        
    return None


if __name__ == '__main__':
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_file = paths.output / f"Table_{chemsys}_Equilibrium_Structures.tex"
    generate_structure_table(paths.data / dbname, chemsys, output_file)

  
