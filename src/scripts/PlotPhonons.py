import paths
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from ase.db import connect
from scipy.interpolate import interp1d

def format_xticklabels(labels, locs):
    # Dictionary to keep track of formatted labels with their locations
    formatted_labels_dict = {}
    
    for loc, label in zip(locs, labels):
        if label == "G":  # Check if label is 'G'
            label = '$\Gamma$'  # Replace with the Greek letter Gamma
        
        # Check if this location already has a label
        if loc in formatted_labels_dict:
            # Append the new label to it with a comma
            formatted_labels_dict[loc] += ', ' + label
        else:
            # Otherwise, just add the label
            formatted_labels_dict[loc] = label
    
    # Now, build the list of formatted labels in the order of locations
    formatted_labels = [formatted_labels_dict[loc] for loc in locs]
    
    return formatted_labels

def normalize_q_points(q_points):
    return (q_points - np.min(q_points)) / (np.max(q_points) - np.min(q_points))


def get_color(value, cmap_name="viridis", vmin=-1.5, vmax=1.5):
    norm = plt.Normalize(vmin=vmin, vmax=vmax)
    cmap = plt.get_cmap(cmap_name)
    return cmap(norm(value))


def PlotPhonons(dbname, model_name, structure_name, num_strains=None):
    """
    This function plots individual phonon bandstructures for different strains and
    then all the strains on one plot with the optical branches semi-transparent.

    TODO: confirm strain is in bars and thus conver to MPa
    """
    with connect(dbname) as db:
        for row in db.select(structure_name=structure_name, model_name=model_name):
            data = row.data
            for strain, pdata in data["strain_phonons"].items():
                fig, ax = plt.subplots()
                for b_idx, energies in pdata["bandstructure"]["band_index"].items():
                    q_points = pdata["bandstructure"]["q_points"]
                    ax.plot(q_points, energies, color="k")

                ax.set_xticks(pdata["bandstructure"]["label_locs"])
                ax.set_xticklabels(format_xticklabels(pdata["bandstructure"]["labels"], pdata["bandstructure"]["label_locs"]))
                _ = [
                    ax.axvline(x=loc, color="gray", linestyle="-", alpha=0.5)
                    for loc in pdata["bandstructure"]["label_locs"]
                ]
                stress = (
                    np.mean(pdata["stress"][:3]) * 0.1
                )  # Confirm units are bar -> MPa
                ax.set_title(f"Strain: {float(strain):.3f}%, Stress: {stress:.3f} MPa")
                ax.set_xlim((0, pdata["bandstructure"]["label_locs"][-1]))
                ax.set_ylabel("Frequency [THz]")
                ax.axhline(y=0.0, color="gray", linestyle="-", alpha=0.5)
                fname = (
                    model_name.upper()
                    + "/"
                    + f"{structure_name}_strain_{float(strain):0.3f}_phonon_bandstructure.png"
                )
                plt.savefig(paths.data / fname)
                plt.close(fig)

    fig, ax = plt.subplots()
    with connect(dbname) as db:
        first_row = db.get(structure_name=structure_name, model_name=model_name)
        f_pdata = next(iter(first_row.data["strain_phonons"].values()))
        labels = f_pdata["bandstructure"]["labels"]
        norm_locs = {}

        #NOTE: I create a evenly space sampling of strains from min to max, but
        # it might also be need to provide ranges based on model
        strains = sorted(first_row.data["strain_phonons"].keys(), key=float)
        if num_strains is not None and num_strains < len(strains):
            strain_indices = np.round(
                np.linspace(0, len(strains) - 1, num_strains)
            ).astype(int)
            selected_strains = [strains[i] for i in strain_indices]
        else:
            selected_strains = strains

        for strain in selected_strains:
            pdata = first_row.data["strain_phonons"][strain]
            q_norm = normalize_q_points(pdata["bandstructure"]["q_points"])
            norm_locs[strain] = normalize_q_points(pdata["bandstructure"]["label_locs"])
            color = get_color(float(strain))

            for b_idx, energies in pdata["bandstructure"]["band_index"].items():
                alpha = 1.0 if b_idx in [0, 1, 2] else 0.375
                label = f"{float(strain):0.2f}%" if b_idx == 0 else ""
                ax.plot(
                    q_norm,
                    energies,
                    color=color,
                    linewidth=1.25,
                    alpha=alpha,
                    label=label,
                )

        f_strain_key = selected_strains[0]
        ax.set_xticks(norm_locs[f_strain_key])
        ax.set_xticklabels(format_xticklabels(labels[: len(norm_locs[f_strain_key])], norm_locs[f_strain_key]))
#        ax.set_xticklabels(labels[: len(norm_locs[f_strain_key])])
        ax.set_xlim((0, norm_locs[f_strain_key][-1]))

        _ = [
            ax.axvline(x=loc, color="gray", linestyle="-", alpha=0.5)
            for loc in norm_locs[f_strain_key]
        ]
        ax.axhline(y=0.0, color="black", linestyle="-", linewidth=0.5, alpha=0.9)
        ax.set_ylabel("Frequency [THz]")

        handles, labels = ax.get_legend_handles_labels()
        labels, handles = zip(
            *sorted(zip(labels, handles), key=lambda t: float(t[0].rstrip("%")))
        )
        ax.legend(
            handles,
            labels,
            loc="upper center",
            bbox_to_anchor=(0.5, -0.10),
            ncol=len(selected_strains),
            fontsize=8,
            fancybox=True,
        )
        plt.subplots_adjust(
            bottom=0.2
        )  # Adjust the bottom to make space for the legend
        plt.savefig(
            paths.figures
            / f"{model_name}_{structure_name}_combined_strain_phonon_bandstructures.png",
            dpi=600,
        )

        return None


# def PlotAllModelPhonons(dbname, models, structure_name, strain=0.00):
#     """
#     This function plots all the models on a single plot at a given strain. The optical branches semi-transparent.

#     """
#     with connect(dbname) as db:
#         fig, ax = plt.subplots()

#         for m, model_name in enumerate(models):
#             first_row = db.get(structure_name=structure_name, model_name=model_name)
#             f_pdata = next(iter(first_row.data["strain_phonons"].values()))
#             labels = f_pdata["bandstructure"]["labels"]
#             norm_locs = {}
#             strains = sorted(first_row.data["strain_phonons"].keys(), key=float)
#             closest_index = np.argmin([abs(float(s) - strain) for s in strains])
#             closest_strain = strains[closest_index]

#             pdata = first_row.data["strain_phonons"][closest_strain]
#             q_norm = normalize_q_points(pdata["bandstructure"]["q_points"])
#             norm_locs = normalize_q_points(pdata["bandstructure"]["label_locs"])
#             color = get_color(m, cmap_name="Dark2",vmin=0,vmax=len(models))

#             for b_idx, energies in pdata["bandstructure"]["band_index"].items():
#                 alpha = 1.0 if b_idx in [0, 1, 2] else 0.375
#                 label = f"{model_name}" if b_idx == 0 else ""
#                 ax.plot(
#                     q_norm,
#                     energies,
#                     color=color,
#                     linewidth=1.25,
#                     alpha=alpha,
#                     label=label,
#                 )

#         ax.set_xticks(norm_locs)
#         ax.set_xticklabels(labels[: len(norm_locs)])
#         ax.set_xlim((0, norm_locs[-1]))

#         _ = [
#             ax.axvline(x=loc, color="gray", linestyle="-", alpha=0.5)
#             for loc in norm_locs
#         ]
#         ax.axhline(y=0.0, color="black", linestyle="-", linewidth=0.5, alpha=0.9)
#         ax.set_ylabel("Frequency [THz]")

#         handles, labels = ax.get_legend_handles_labels()
#         #labels, handles = zip(
#         #    *sorted(zip(labels, handles), key=lambda t: float(t[0].rstrip("%")))
#         #)
#         ax.legend(
#             handles,
#             labels,
#             loc="upper center",
#             bbox_to_anchor=(0.5, -0.10),
#             ncol=len(models),
#             fontsize=8,
#             fancybox=True,
#         )
#         plt.subplots_adjust(
#             bottom=0.2
#         )  # Adjust the bottom to make space for the legend
#         plt.savefig(
#             paths.figures
#             / f"{structure_name}_combined_models_phonon_bandstructures.png",
#             dpi=600,
#         )

#         return None



def PlotAllModelPhonons(dbname, models, structure_name, strain=0.00, columns=3,fade_optical=False):
    """
    This function creates a grid of phonon bandstructure plots for each model.
    """
    with connect(dbname) as db:
        # Determine the grid size
        rows = int(np.ceil(len(models) / columns))
        
        # Create a grid of subplots
        fig = plt.figure(figsize=(12, 4 * rows))
        gs = GridSpec(rows, columns, figure=fig)

        for m, model_name in enumerate(models):
            ax = fig.add_subplot(gs[m // columns, m % columns])

            first_row = db.get(structure_name=structure_name, model_name=model_name)
            #NOTE: we skip B19 and B19P for ALIGNN due to performance issues
            #if model_name == "ALIGNN":
            #    if structure_name != "B2":
            #        UserWarning(f"Skipping: {structure_name} {model_name}")
            #        continue
                
            f_pdata = next(iter(first_row.data["strain_phonons"].values()))
            labels = f_pdata["bandstructure"]["labels"]
            strains = sorted(first_row.data["strain_phonons"].keys(), key=float)
            closest_index = np.argmin([abs(float(s) - strain) for s in strains])
            closest_strain = strains[closest_index]

            pdata = first_row.data["strain_phonons"][closest_strain]
            q_norm = normalize_q_points(pdata["bandstructure"]["q_points"])
            norm_locs = normalize_q_points(pdata["bandstructure"]["label_locs"])

            for b_idx, energies in pdata["bandstructure"]["band_index"].items():
                if fade_optical:
                    alpha = 1.0 if b_idx in [0, 1, 2] else 0.375
                else:
                    alpha = 1.0
                ax.plot(
                    q_norm,
                    energies,
                    color='black',
                    linewidth=1.25,
                    alpha=alpha
                )

            ax.set_xticks(norm_locs)
            ax.set_xticklabels(format_xticklabels(labels[: len(norm_locs)], norm_locs))
            #ax.set_xticklabels(labels[: len(norm_locs)])
            ax.set_xlim((0, norm_locs[-1]))
            ax.axhline(y=0.0, color="black", linestyle="-", linewidth=0.5, alpha=0.9)
            ax.set_ylabel("Frequency [THz]")
            ax.set_title(model_name, fontsize=10)

            # Add subplot label
            ax.text(-0.05, 1.05, f"{chr(97 + m)})", transform=ax.transAxes, size=12, weight='bold')

            _ = [
                ax.axvline(x=loc, color="gray", linestyle="-", alpha=0.5)
                for loc in norm_locs
            ]
        
        # Adjust layout
        plt.tight_layout()
        plt.savefig(
            paths.figures
            / f"{structure_name}_combined_models_phonon_bandstructures.png",
            dpi=600,
        )

        return None
