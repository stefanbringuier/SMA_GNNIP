import sys
import paths
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
import argparse

from ase.db import connect


def format_xticklabels(labels, locs):
    # Dictionary to keep track of formatted labels with their locations
    formatted_labels_dict = {}

    for loc, label in zip(locs, labels):
        if label == "G":  # Check if label is 'G'
            label = "$\Gamma$"  # Replace with the Greek letter Gamma

        # Check if this location already has a label
        if loc in formatted_labels_dict:
            # Append the new label to it with a comma
            formatted_labels_dict[loc] += ", " + label
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


def plot_default_phonons(bs, dos, info, ymaxlim=9.0, of="."):
    """
    Plot results from ASE phonon calculation.
    """
    fig = plt.figure(1, figsize=(7, 4))
    ax = fig.add_axes([0.12, 0.07, 0.67, 0.85])
    dosax = fig.add_axes([0.8, 0.07, 0.17, 0.85])

    xticks = bs.get_labels()[1].tolist()
    xlabels = bs.get_labels()[2]

    for xt in xticks:
        ax.axvline(x=xt, color="lightgray")
    ax.axhline(y=0.0, color="black")

    q_points = bs.get_labels()[0]
    energies_THz = EV_to_THz * bs.energies
    for n in range(len(energies_THz[0, 1, :])):
        omega = energies_THz[0, :, n]
        ax.plot(q_points, omega, "k-", lw=2)

    ax.set_title(
        f"Volume: {info['vol']:.3f}Å^3\
        Strain: {info['strain']:.2f}% \
        Potential: {info['potname']}",
        fontsize=10,
    )

    ax.set_xticks(xticks)
    ax.set_xticklabels(xlabels)

    ax.set_xlim((0.0, q_points[-1]))
    ax.set_ylim((-2.0, ymaxlim))
    ax.set_ylabel("Frequency ($\mathrm{THz}$)", fontsize=22)

    dosax.fill_between(
        dos.get_weights(),
        EV_to_THz * dos.get_energies(),
        y2=0,
        color="grey",
        edgecolor="k",
        lw=1,
    )

    ymaxlim = np.max(energies_THz) * 1.05
    dosax.set_ylim((-2.0, ymaxlim))
    dosax.set_yticks([])
    dosax.set_xticks([])
    dosax.set_xlabel("DOS", fontsize=18)
    fig.savefig(
        os.path.join(
            of,
            f"{info['chemsys']}_{info['structure']}_{info['potname']}_{info['strain']:.3f}_Phonons.png",
        )
    )
    fig.clear()

    return None


def plot_individual_phonon_bandstructure(dbname, chemsys, model_name, structure_name):
    """
    Plot individual phonon bandstructures for different strains from ASE database.

    Args:
        dbname (str): Name of the database.
        chemsys (str): Chemical system identifier.
        model_name (str): Name of the model.
        structure_name (str): Name of the structure.

    Returns:
        None
    """
    with connect(dbname) as db:
        for row in db.select(structure_name=structure_name, model_name=model_name):
            if chemsys != "".join(sorted(set(row.symbols))):
                break

            data = row.data
            for strain, pdata in data["strain_phonons"].items():
                fig, ax = plt.subplots()
                for b_idx, energies in pdata["bandstructure"]["band_index"].items():
                    q_points = pdata["bandstructure"]["q_points"]
                    ax.plot(q_points, energies, color="k")

                ax.set_xticks(pdata["bandstructure"]["label_locs"])
                ax.set_xticklabels(
                    format_xticklabels(
                        pdata["bandstructure"]["labels"],
                        pdata["bandstructure"]["label_locs"],
                    )
                )
                # Calculate midpoints for q_vec_labels
                midpoints = 0.5 * np.array(
                    pdata["bandstructure"]["label_locs"][:-1]
                    + pdata["bandstructure"]["label_locs"][1:]
                )
                for midpoint, vec_label in zip(
                    midpoints, pdata["bandstructure"]["q_vec_labels"]
                ):
                    ax.text(
                        midpoint,
                        ax.get_ylim()[1],
                        vec_label,
                        ha="center",
                        va="bottom",
                        transform=ax.transData,
                        size=10,
                    )

                _ = [
                    ax.axvline(x=loc, color="gray", linestyle="-", alpha=0.5)
                    for loc in pdata["bandstructure"]["label_locs"]
                ]
                stress = (
                    np.mean(pdata["stress"][:3]) * 160.21766208  # eV/Ang^3 to GPa
                )  # WARNING: Need to confirm units are eV/Ang^3 for all calculators
                ax.set_title(f"Strain: {float(strain):.3f}%, Stress: {stress:.3f} GPa")
                ax.set_xlim((0, pdata["bandstructure"]["label_locs"][-1]))
                ax.set_ylabel("Frequency [THz]")
                ax.axhline(y=0.0, color="gray", linestyle="-", alpha=0.5)
                fname = (
                    model_name.upper()
                    + "/"
                    + f"{chemsys}_{structure_name}_strain_{float(strain):0.3f}_phonon_bandstructure.png"
                )
                plt.savefig(paths.data / fname)
                plt.close(fig)

    return None


def plot_all_strains_phonons(
    dbname, chemsys, model_name, structure_name, num_strains=None
):
    """
    Plot all strain phonon bandstructures on one plot with optical branches semi-transparent.

    Args:
        dbname (str): Name of the database.
        chemsys (str): Chemical system identifier.
        model_name (str): Name of the model.
        structure_name (str): Name of the structure.
        num_strains (int, optional): Number of strains to plot. Defaults to None.

    Returns:
        None
    """
    fig, ax = plt.subplots()

    with connect(dbname) as db:
        first_row = db.get(structure_name=structure_name, model_name=model_name)
        if chemsys != "".join(sorted(set(first_row.symbols))):
            return None

        f_pdata = next(iter(first_row.data["strain_phonons"].values()))
        labels = f_pdata["bandstructure"]["labels"]
        norm_locs = {}

        # NOTE: I create a evenly space sampling of strains from min to max, but
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
        ax.set_xticklabels(
            format_xticklabels(
                labels[: len(norm_locs[f_strain_key])], norm_locs[f_strain_key]
            )
        )

        # Calculate midpoints for q_vec_labels
        midpoints = 0.5 * (norm_locs[f_strain_key][:-1] + norm_locs[f_strain_key][1:])
        q_vec_labels = f_pdata["bandstructure"]["q_vec_labels"]
        for midpoint, vec_label in zip(midpoints, q_vec_labels):
            ax.text(
                midpoint,
                ax.get_ylim()[1],
                vec_label,
                ha="center",
                va="bottom",
                transform=ax.transData,
                size=6,
            )

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
            / f"{chemsys}_{model_name}_{structure_name}_StrainsPhononBandstructures.png",
            dpi=600,
        )

    return None


def plot_all_model_phonons(
    dbname, chemsys, models, structure_name, strain=0.00, columns=3, fade_optical=False
):
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

            if chemsys != "".join(sorted(set(first_row.symbols))):
                break

            f_pdata = next(iter(first_row.data["strain_phonons"].values()))
            labels = f_pdata["bandstructure"]["labels"]
            strains = sorted(first_row.data["strain_phonons"].keys(), key=float)
            closest_index = np.argmin([abs(float(s) - strain) for s in strains])
            closest_strain = strains[closest_index]

            pdata = first_row.data["strain_phonons"][closest_strain]
            q_norm = normalize_q_points(pdata["bandstructure"]["q_points"])
            norm_locs = normalize_q_points(pdata["bandstructure"]["label_locs"])

            # NOTE: Used for B2
            for b_idx, energies in pdata["bandstructure"]["band_index"].items():
                if fade_optical:
                    alpha = 1.0 if b_idx in [0, 1, 2] else 0.375
                else:
                    alpha = 1.0
                ax.plot(q_norm, energies, color="black", linewidth=1.25, alpha=alpha)

            ax.set_xticks(norm_locs)
            ax.set_xticklabels(format_xticklabels(labels[: len(norm_locs)], norm_locs))
            # Calculate midpoints for q_vec_labels
            midpoints = 0.5 * (norm_locs[:-1] + norm_locs[1:])
            for midpoint, vec_label in zip(
                midpoints, pdata["bandstructure"]["q_vec_labels"]
            ):
                ax.text(
                    midpoint,
                    ax.get_ylim()[1],
                    vec_label,
                    ha="center",
                    va="bottom",
                    transform=ax.transData,
                    size=6,
                )

            # ax.set_xticklabels(labels[: len(norm_locs)])
            ax.set_xlim((0, norm_locs[-1]))
            ax.axhline(y=0.0, color="black", linestyle="-", linewidth=0.5, alpha=0.9)
            ax.set_ylabel("Frequency [THz]")
            # ax.set_title(model_name, fontsize=10)

            # Add subplot label
            ax.text(
                -0.05,
                1.05,
                f"{chr(97 + m)})",
                transform=ax.transAxes,
                size=12,
                weight="bold",
            )

            _ = [
                ax.axvline(x=loc, color="gray", linestyle="-", alpha=0.5)
                for loc in norm_locs
            ]

        # Adjust layout
        plt.tight_layout()
        plt.savefig(
            paths.figures
            / f"{chemsys}_{structure_name}_ModelsPhononBandstructures.png",
            dpi=600,
        )

        return None

def parse_arguments():
    parser = argparse.ArgumentParser(description="Process structures and models.")
    parser.add_argument("dbname", type=str, help="Database name")
    parser.add_argument("num_strains", type=int, help="Number of strains")
    parser.add_argument("chemsys", type=str, help="Chemical system")
    parser.add_argument("--structures", nargs='+', type=str, required=True, help="List of structures")
    parser.add_argument("--models", nargs='+', type=str, required=True, help="List of models")
    return parser.parse_args()

if __name__ == "__main__":
    args = parse_arguments()
    dbname = paths.data / args.dbname 
    for s in args.structures:
        plot_all_model_phonons(dbname, args.chemsys, args.models, s)
    for m in args.models:
        for s in args.structures:
            plot_all_strains_phonons(dbname, args.chemsys, m, s, num_strains=args.num_strains)
