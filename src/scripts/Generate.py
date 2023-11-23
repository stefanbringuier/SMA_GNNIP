import warnings
import argparse
import os
import paths

from NiTiStructures import *
from Calculators import *

from CalculateEOS import *
from CalculateElastic import *
from CalculatePhonons import *

# Surpress warning from dependencies for which I can do nothing about.
warnings.filterwarnings('ignore', category=UserWarning, module='ase.spacegroup')
warnings.filterwarnings("ignore", category=UserWarning, module="dgl.backend.pytorch.tensor")

def main():
    '''
     Provides equation of state calculations of NiTi. Be mindful that correct python environment is
     activated for selected model (i.e., ASE Calculator)
     '''
    parser = argparse.ArgumentParser(description="""Equation of State Calculations of NiTi.
    Be mindful that the correct python environment is activated for the model (i.e., ASE Calculator) you
    intend to use!!""")
    parser.add_argument('--dbname', type=str, default='NiTi_EOS.json', help='ASE Database file name')
    parser.add_argument('--dbfolder', type=str, default='.', help='ASE Database folder')
    parser.add_argument('--model', choices=['Mutter', 'Zhong', 'M3GNet', 'CHGNet','MACE','ALIGNN'], default='Mutter', help='Choose ASECalculator Model')
    parser.add_argument('--structure',default='B2')
    parser.add_argument('--n',type=int,default=40,help='Number of configurations in strain scan of EOS')
    parser.add_argument('--nph',type=int,default=13,help='Number of configurations in strain scan of Phonons')
    parser.add_argument('--min_strain', type=float, default=-0.11, help='Minimum strain, fractional')
    parser.add_argument('--max_strain', type=float, default=0.11, help='Maximum strain, fractional')
    parser.add_argument('--eos_fit', type=str, default='sj', help='EOS fit function')
    args = parser.parse_args()
    

    outfolder = paths.data / args.model.upper()
    os.makedirs(outfolder, exist_ok=True)

    asecalc = GetCalculator(args)
    structure = GetStructure(args.structure)
    structure.calc = asecalc
    structure.info["model_name"] = args.model
    
    
    db_path_file = paths.data / args.dbname
    print(f"""
    -----------------------------------------------------
    EOS Calculation of {args.structure} using {args.model}
    -----------------------------------------------------
    """)
    CalculateEOS(structure,
             (args.model, asecalc),
             outfolder,
             strain=np.linspace(args.min_strain,args.max_strain,args.n),
             dbname=db_path_file,
             eos_eq=args.eos_fit)

    print(f"""
    -----------------------------------------------------
    Phonon Calculation of {args.structure} using {args.model}
    -----------------------------------------------------
    """)
    strain = np.linspace(args.min_strain,args.max_strain,args.nph)
    CalculatePhonons(db_path_file,args.structure,(args.model,asecalc),strain=strain)
    
    # NOT VALIDATED
    #CalculateElasticConstants(db_path_file,(args.model,asecalc))
    
    
if __name__ == "__main__":
    main()
