import sys

import numpy as np
import os
import paths


import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator, AutoLocator
from cycler import cycler
import seaborn as sns  # For color palette
from ase.db import connect

SUBPLOT_ORDER = {'a': 'Mutter',
                 'b': 'Zhong',
                 'c': 'M3GNet',
                 'd': 'CHGNet',
                 'e': 'MACE',
                 'f': 'ALIGNN',
                 }

def plot_eos(db_path, figsize=(3.4, 5.75), dpi=600, xlim=(0.75, 1.25), ignore_list=['B32', 'R_Phase'], 
             color_palette="muted", n_colors=10, wspace=0.3, hspace=0.3, linewidth=1, linestyle='-', 
             fontsize=8, tick_fontsize=6, legend_fontsize=8, legend_ncol=3):
   
    # Connect to the ASE database
    db = connect(db_path)
    
    # Gather all unique calculator models and structure names
#    calculator_models = set(row.model_name for row in db.select())

    calculator_models = sorted(set(row.model_name for row in db.select()),
                               key=lambda model: list(SUBPLOT_ORDER.values()).index(model))

    structure_names = set(row.structure_name for row in db.select())

    
    # Create a figure with a 2x3 grid of subplots
    fig, axs = plt.subplots(3, 2, figsize=figsize)
    axs = axs.ravel()  # Flatten the array of axes for easier indexing

    # Create a color cycle iterator using a slightly modified color palette
    color_cycle = sns.color_palette(color_palette, n_colors=n_colors)
    
    label_file_path = paths.figures / 'NiTi_EOS_Comparison_label.txt'
    with open(label_file_path, 'w') as label_file:
        # Iterate over calculator models
        for i, calculator_model in enumerate(calculator_models):

            if calculator_model not in SUBPLOT_ORDER.values():
                raise ValueError(f"{calculator_model} model is not expected")
            
            ax = axs[i]  # Select the current axis
            alphabetic_label = f'({chr(97 + i)})'
            ax.set_title(alphabetic_label, loc='left', fontsize=fontsize)  # Set title with alphabetic label
            label_file.write(f'{alphabetic_label}: {calculator_model}\n')  # Write label to file
            
            # Set the color cycle for the current axis
            ax.set_prop_cycle(cycler(color=color_cycle))
            
            ax.set_xlim(xlim)
            #ax.text(-5.0, 1.1, f'{calculator_model}', ha='center', fontsize=tick_fontsize)  # Single x-label
            #ax.xaxis.set_major_locator(ticker.MaxNLocator(nbins=6))
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

            ax.tick_params(axis='both',direction='in',length=2)
            ax.tick_params(labelsize=tick_fontsize)

            ax.axvline(x=1.0, color='lightgray', linestyle=':', linewidth=linewidth)
            
            # Iterate over structure names
            for structure_name in structure_names:
                if structure_name in ignore_list:
                    continue  # Skip the structures in the ignore list
                
                for row in db.select(model_name=calculator_model, structure_name=structure_name):
                    volumes = np.array(row.data['volumes'])
                    energies = np.array(row.data['energies'])
                    #NOTE: We can't use EOS fit volume, as the fits aren't always that good.
                    #vi = row.v0
                    min_energy_index = np.argmin(energies)
                    vi = volumes[min_energy_index]
                    natoms = row.natoms
                    
                    vv_ratio = volumes / vi
                    energy_per_atom = energies / natoms
                    ax.plot(vv_ratio, energy_per_atom, linewidth=linewidth, linestyle=linestyle, label=structure_name)

                    # Update the global y-axis limits
                    ymin =  min(energy_per_atom) - 0.1
                    ymax =  max(energy_per_atom) * 1.05
                    ax.set_ylim((ymin, ymax))

        handles, labels = axs[-1].get_legend_handles_labels()
        fig.legend(handles, labels, loc='upper center', ncol=legend_ncol, 
                   fontsize=legend_fontsize, edgecolor='black', 
                   fancybox=True)

        fig.text(0.5, 0.05, 'V/V$_o$', ha='center', fontsize=fontsize)  # Single x-label
        fig.text(0.01, 0.5, 'Energy per atom (eV)', va='center', rotation='vertical', fontsize=fontsize)  # Single y-label
        plt.subplots_adjust(wspace=wspace, hspace=hspace)  # Adjust spacing as needed

        file_name = 'NiTi_EOS_Comparison.png'
        plt.savefig(paths.figures / file_name, dpi=dpi)  # Increased resolution
        print(f'Saved plot to {os.path.abspath(file_name)}')
        print(f'Saved labels to {os.path.abspath(label_file_path)}')  # Debug statement

if __name__ == "__main__":
    db_path = paths.data / sys.argv[1]  # Specify the path to your ASE database
    plot_eos(db_path)
