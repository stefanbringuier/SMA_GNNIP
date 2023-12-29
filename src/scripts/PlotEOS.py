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

from Config import SUBPLOT_ORDER

def plot_eos(
    dbname,
    output_plot,
    figsize=(3.4, 5.75),
    dpi=600,
    xlim=(0.75, 1.25),
    ignore_list=["B32", "R_Phase"],
    color_palette="muted",
    n_colors=10,
    wspace=0.3,
    hspace=0.3,
    linewidth=1,
    linestyle="-",
    fontsize=8,
    tick_fontsize=6,
    legend_fontsize=8,
    legend_ncol=5,
    label_notes="EOS_Labels.txt",
):

    db = connect(dbname)

    calculator_models = sorted(
        set(row.model_name for row in db.select()),
        key=lambda model: list(SUBPLOT_ORDER.values()).index(model),
    )

    structure_names = set(row.structure_name for row in db.select())

    # Create a figure with a 2x3 grid of subplots
    fig, axs = plt.subplots(3, 2, figsize=figsize)
    axs = axs.ravel()  # Flatten the array of axes for easier indexing

    # Create a color cycle iterator using a slightly modified color palette
    color_cycle = sns.color_palette(color_palette, n_colors=n_colors)

    label_file_path = paths.figures / label_notes
    with open(label_file_path, "w") as label_file:
        for i, calculator_model in enumerate(calculator_models):
            if calculator_model not in SUBPLOT_ORDER.values():
                raise ValueError(f"{calculator_model} model is not expected")

            ax = axs[i] 
            alphabetic_label = f"({chr(97 + i)})"
            ax.set_title(
                alphabetic_label, loc="left", fontsize=fontsize
            ) 
            label_file.write(
                f"{alphabetic_label}: {calculator_model}\n"
            )

          
            ax.set_prop_cycle(cycler(color=color_cycle))
            ax.set_xlim(xlim)
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=6))

            ax.tick_params(axis="both", direction="in", length=2)
            ax.tick_params(labelsize=tick_fontsize)

            ax.axvline(x=1.0, color="lightgray", linestyle=":", linewidth=linewidth)

            for structure_name in structure_names:
                if structure_name in ignore_list:
                    continue  # Skip the structures in the ignore list

                for row in db.select(
                    model_name=calculator_model, structure_name=structure_name
                ):
                    volumes = np.array(row.data["volumes"])
                    energies = np.array(row.data["energies"])
                    # NOTE: We can't use EOS fit volume, as the fits aren't always that good.
                    # vi = row.v0
                    min_energy_index = np.argmin(energies)
                    vi = volumes[min_energy_index]
                    natoms = row.natoms

                    vv_ratio = volumes / vi
                    energy_per_atom = energies / natoms
                    ax.plot(
                        vv_ratio,
                        energy_per_atom,
                        linewidth=linewidth,
                        linestyle=linestyle,
                        label=structure_name,
                    )

                    # Update the global y-axis limits
                    ymin = min(energy_per_atom) - 0.1
                    ymax = max(energy_per_atom) * 1.05
                    ax.set_ylim((ymin, ymax))

        handles, labels = axs[-1].get_legend_handles_labels()
        fig.legend(
            handles,
            labels,
            loc="upper center",
            ncol=legend_ncol,
            fontsize=legend_fontsize,
            edgecolor="black",
            fancybox=True,
        )

        # Single x-axis label
        fig.text(0.5, 0.05, "V/V$_o$", ha="center", fontsize=fontsize)
        # Single y-axis label
        fig.text(
            0.01,
            0.5,
            "Energy per atom (eV)",
            va="center",
            rotation="vertical",
            fontsize=fontsize,
        ) 
        plt.subplots_adjust(wspace=wspace, hspace=hspace)
        plt.savefig(output_plot, dpi=dpi)

        return None


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_plot = paths.figures / f"{chemsys}_EquationOfStates.png"
    plot_eos(paths.data / dbname, paths.figures / output_plot)
