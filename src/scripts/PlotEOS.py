import sys
import warnings

import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np
import paths
import scienceplots
import seaborn as sns  # For color palette
from ase.db import connect
from cycler import cycler
from matplotlib.ticker import MultipleLocator

plt.style.use(["science", "scatter", "no-latex"])

from PlotConfigs import SUBPLOT_ORDER


def plot_eos(
    dbname,
    chemsys,
    output_plot,
    figsize=(7.24, 7.24),
    dpi=600,
    xlim=(0.75, 1.25),
    ignore_list=["B32", "R_Phase"],
    color_palette="muted",
    n_colors=10,
    wspace=0.3,
    hspace=0.3,
    linewidth=0.75,
    linestyle="-",
    fontsize=5,
    tick_fontsize=5,
    legend_fontsize=6,
    legend_ncol=5,
    label_notes="EOS_Labels.txt",
):
    db = connect(dbname)

    # Select unique models if in SUBPLOT_ORDER and sort
    plot_models = SUBPLOT_ORDER[chemsys].values()

    calculator_models = sorted(
        set(
            row.model_name
            for row in db.select(chemsys=chemsys)
            if row.model_name in plot_models
        ),
        key=lambda model: list(SUBPLOT_ORDER[chemsys].values()).index(model),
    )
    structure_names = set(row.structure_name for row in db.select(chemsys=chemsys))

    # Create a figure with a 2x3 grid of subplots
    fig, axs = plt.subplots(2, 3)  # , figsize=figsize)
    axs = axs.ravel()  # Flatten the array of axes for easier indexing

    # Create a color cycle iterator using a slightly modified color palette
    color_cycle = sns.color_palette(color_palette, n_colors=n_colors)

    label_file_path = paths.figures / label_notes
    with open(label_file_path, "w") as label_file:
        for i, calculator_model in enumerate(calculator_models):
            if calculator_model not in SUBPLOT_ORDER[chemsys].values():
                warnings.warn(
                    f"{calculator_model} model results not used/expected", UserWarning
                )

            ax = axs[i]
            alphabetic_label = f"({chr(97 + i)})"
            ax.set_title(alphabetic_label, loc="left", fontsize=fontsize)
            label_file.write(f"{alphabetic_label}: {calculator_model}\n")

            ax.set_prop_cycle(cycler(color=color_cycle))
            ax.set_xlim(xlim)
            ax.xaxis.set_major_locator(MultipleLocator(0.1))
            ax.yaxis.set_major_locator(ticker.MaxNLocator(nbins=4))

            # ax.tick_params(axis="both", direction="in", length=2)
            # ax.tick_params(labelsize=tick_fontsize)

            ax.axvline(x=1.0, color="lightgray", linestyle=":", linewidth=linewidth)

            for structure_name in structure_names:
                if structure_name in ignore_list:
                    continue  # Skip the structures in the ignore list

                ymin, ymax = [], []
                for row in db.select(
                    chemsys=chemsys,
                    model_name=calculator_model,
                    structure_name=structure_name,
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
                    ymin.append(min(energy_per_atom))
                    ymax.append(max(energy_per_atom))

                ax.set_ylim((min(ymin) - 0.05, max(ymax) * 1.10))
                ax.tick_params(axis="x", which="minor", bottom=False, top=False)
                # ax.tick_params(axis='y', which='minor', right=False, left=False)
                ax.tick_params(axis="x", labelsize=fontsize)
                ax.tick_params(axis="y", labelsize=fontsize)

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
        plt.close()
        return None


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_plot = paths.figures / f"{chemsys}_EquationOfStates.png"
    plot_eos(paths.data / dbname, chemsys, paths.figures / output_plot)
