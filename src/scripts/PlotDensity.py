import sys

import matplotlib.pyplot as plt
import paths
import scienceplots
from matplotlib.ticker import MaxNLocator

plt.style.use(["science", "scatter", "no-latex"])

from ase.db import connect
from PlotConfigs import ORDER, STRUCTURE_ORDER


def plot_density(dbname, chemsys, output_filename):
    db = connect(dbname)

    sorder = STRUCTURE_ORDER[chemsys]
    default_order = max(sorder.values())
    unique_structure_types = set(entry.structure_name for entry in db.select())
    structure_order = sorted(
        unique_structure_types, key=lambda x: sorder.get(x, default_order)
    )

    # Data containers
    structure_names = []
    model_names = set()
    density_data = {}

    # Collecting data
    for name in structure_order:
        entries = db.select(chemsys=chemsys, structure_name=name)
        structure_names.append(name.replace("_", "-"))
        density_data[name] = {}

        for entry in entries:
            model_name = entry.get("model_name", "NA")
            density = entry.get("density", 0)  # Default to 0 if not available
            if model_name in ORDER[chemsys]:
                density_data[name][model_name] = density
                model_names.add(model_name)

    # structure_names = sorted(structure_names)
    model_names = sorted(model_names, key=lambda x: ORDER[chemsys].get(x, float("inf")))

    # Fancy colors and patterns
    colors = ["red", "green", "blue", "cyan", "magenta", "yellow", "black"]
    patterns = ["//", "+", "*", "x", "o", "-", "O", "||", "++"]

    # Plotting
    fig, ax = plt.subplots()
    bar_width = 0.125
    space_between_bars = 0.075  # Adjust this value as needed
    indices = list(range(len(structure_names)))
    num_models = len(model_names)
    group_centers = [
        1.65 * x + (num_models / 2) * (bar_width + space_between_bars) for x in indices
    ]

    for i, model_name in enumerate(model_names):
        offset = [(1.5 * x + i * (bar_width + space_between_bars)) for x in indices]
        density_values = [
            density_data[struct].get(model_name, 0) for struct in structure_order
        ]
        # Assign color and pattern
        color = colors[i % len(colors)]
        pattern = patterns[i % len(patterns)]
        ax.bar(
            offset,
            density_values,
            bar_width,
            label=model_name,
            color=color,
            edgecolor="black",
            hatch=pattern,
            alpha=0.75,
        )

    ax.set_xlabel("Structure Name")
    ax.set_ylabel("Density [g/cm$^3$]")
    ax.set_xticks(group_centers)
    xlabels = [x.replace("B19P", "B19'") for x in structure_names]
    ax.set_xticklabels(xlabels, rotation=0, ha="right")
    ax.tick_params(axis="x", which="both", length=0)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim((5, None))
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=6)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_plot = paths.figures / f"{chemsys}_Equilibrium_Density.png"
    plot_density(paths.data / sys.argv[1], chemsys, output_plot)
