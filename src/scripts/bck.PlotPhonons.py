import paths
import numpy as np
import matplotlib.pyplot as plt
from ase.db import connect


def normalize_q_points(q_points):
    return (q_points - np.min(q_points)) / (np.max(q_points) - np.min(q_points))

def get_color(value, cmap_name='viridis', vmin=-0.12, vmax=0.03):
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)
    return cmap(norm(value))
#def get_color(value, cmap_name='Blues', vmin=-0.11, vmax=0.01, start=0.1, stop=0.9):
#    norm = plt.Normalize(vmin=vmin, vmax=vmax)
#    cmap = plt.get_cmap(cmap_name)
#    scaled_value = start + (stop - start) * norm(value)
#    return cmap(scaled_value)

def PlotPhonons(dbname, model_name, structure_name):
    with connect(dbname) as db:
        for row in db.select(structure_name=structure_name, model_name=model_name):
            data = row.data
            for strain, pdata in data['strain_phonons'].items():
                fig, ax = plt.subplots()
                for b_idx, energies in pdata['bandstructure']['band_index'].items():
                    q_points = pdata['bandstructure']['q_points']
                    ax.plot(q_points, energies, color='k')
                
                ax.set_xticks(pdata['bandstructure']['label_locs'])
                ax.set_xticklabels(pdata['bandstructure']['labels'])
                _ = [ax.axvline(x=loc, color='gray', linestyle='-', alpha=0.5) for loc in pdata['bandstructure']['label_locs']]
                stress = np.mean(pdata["stress"][:3]) * 0.1  # Confirm units are bar -> MPa
                ax.set_title(f'Strain: {float(strain)*100}%, Stress: {stress:.3f} MPa')
                ax.set_xlim((0, pdata['bandstructure']['label_locs'][-1]))
                ax.set_ylabel("Frequency [THz]")
                plt.savefig(f'/tmp/{structure_name}_phonon_bandstructure_{strain}.png')
                plt.close(fig)
    
    fig, ax = plt.subplots()
    with connect(dbname) as db:
        first_row = db.get(structure_name=structure_name, model_name=model_name)
        f_pdata = next(iter(first_row.data['strain_phonons'].values()))
        labels = f_pdata['bandstructure']['labels']
        norm_locs = {}

        for strain, pdata in first_row.data['strain_phonons'].items():
            q_norm = normalize_q_points(pdata['bandstructure']['q_points'])
            norm_locs[strain] = normalize_q_points(pdata['bandstructure']['label_locs'])
            color = get_color(float(strain))

            for b_idx, energies in pdata['bandstructure']['band_index'].items():
                alpha = 1.0 if b_idx in [0, 1, 2] else 0.2
                label = f"{float(strain)*100:0.2f}%" if b_idx == 0 else ""
                ax.plot(q_norm, energies, color=color, linewidth=1.25, alpha=alpha, label=label)
    
        f_strain_key = list(first_row.data['strain_phonons'].keys())[0]
        ax.set_xticks(norm_locs[f_strain_key])
        ax.set_xticklabels(labels[:len(norm_locs[f_strain_key])])
        ax.set_xlim((0, norm_locs[f_strain_key][-1]))

        _ = [ax.axvline(x=loc, color='gray', linestyle='-', alpha=0.5) for loc in norm_locs[f_strain_key]]
        ax.axhline(y=0.0, color='black', linestyle='-', linewidth=0.5, alpha=0.9)
        
        ax.set_ylabel("Frequency [THz]")
        #Sort legend
        handles, labels = ax.get_legend_handles_labels()
        labels, handles = zip(*sorted(zip(labels, handles), key=lambda t: float(t[0].rstrip('%'))))
        legend = ax.legend(handles,labels,loc='upper center', bbox_to_anchor=(0.5, -0.10), ncol=9, fontsize=8, fancybox=True)
        plt.subplots_adjust(bottom=0.2)  # Adjust the bottom to make space for the legend
        plt.savefig(f'/tmp/{model_name}_{structure_name}_combined_phonon_bandstructures.png', dpi=300)

# Usage example:
PlotPhonons(paths.data / 'NiTi_Structures.json', 'Mutter', 'B2')
