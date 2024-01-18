import sys

import matplotlib.pyplot as plt
import paths
import scienceplots
from matplotlib.ticker import MaxNLocator

plt.style.use(["science", "scatter", "no-latex"])

from ase.db import connect
from PlotConfigs import ORDER, STRUCTURE_ORDER


def get_exp_B2_ecoh(chemsys):
    """
    Data point for NiTi digitally extracted from Fig. 2 in Vandermause, J.; et al.  http://arxiv.org/abs/2401.05568.
    """

    ecoh = {
        "NiTi": (-5.0125, "Exp."),
        "PtTi": (None, "Exp."),
    }
    return ecoh[chemsys]


def get_materials_project_B2_ecoh(chemsys):
    """Results coded from Materials Project.

    The obtained valued for B2 structure from  https://materialsproject.com/materials/mp-571

    Code is provided in this doc string as the mp_api pip package is required along with a api key.

    ```python
    with MPRester(mp_apikey) as mpr:
         docs = mpr.summary.search(material_ids=["mp-571"])
    docs[0].energy_per_atom
    ```

    Provides the corrected energy per atom in eV.

    Notes:
    - The `energy_per_atom` for PtTi appears wrong (-33.08 eV) use `uncorrected_energy_per_atom`
    """
    ecoh = {
        "NiTi": (-7.185175e0, "mp-571"),
        "PtTi": (-7.767860005e0, "mp-11552"),
    }

    return ecoh[chemsys]


def plot_cohesive_energy(dbname, chemsys, output_filename):
    db = connect(dbname)

    sorder = STRUCTURE_ORDER[chemsys]
    default_order = max(sorder.values())
    unique_structure_types = set(entry.structure_name for entry in db.select())
    structure_order = sorted(
        unique_structure_types, key=lambda x: sorder.get(x, default_order)
    )

    # unique_structure_types = set(
    #     entry.structure_name for entry in db.select(chemsys=chemsys)
    # )

    # Data containers
    structure_names = []
    model_names = set()
    ecoh_data = {}

    # Collecting data
    for name in structure_order:
        entries = db.select(chemsys=chemsys, structure_name=name)
        structure_names.append(name.replace("_", "-"))
        ecoh_data[name] = {}

        for entry in entries:
            model_name = entry.get("model_name", "NA")
            ecoh = entry.get("ecoh", 0)  # Default to 0 if not available
            if model_name in ORDER[chemsys]:
                ecoh_data[name][model_name] = ecoh
                model_names.add(model_name)

    # structure_names = sorted(structure_names)
    model_names = sorted(model_names, key=lambda x: ORDER[chemsys].get(x, float("inf")))

    # Fancy colors and patterns
    colors = ["red", "green", "blue", "cyan", "magenta", "yellow", "black"]
    patterns = ["//", "+", "*", "x", "o", "||", "O", "-", "\\"]

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
        # offset = [(x + i * bar_width) for x in indices]
        offset = [(1.5 * x + i * (bar_width + space_between_bars)) for x in indices]
        ecoh_values = [
            ecoh_data[struct].get(model_name, 0) for struct in structure_order
        ]
        # Assign color and pattern
        color = colors[i % len(colors)]
        pattern = patterns[i % len(patterns)]
        ax.bar(
            offset,
            ecoh_values,
            bar_width,
            label=model_name,
            color=color,
            edgecolor="black",
            hatch=pattern,
            alpha=0.75,
        )

    # TODO: Need to manually add ab-initio cohesive energy results.
    # Most papers are shoring difference B2-Structure in meV or
    # formation energies. For now I will use the Materials Project
    # B2 structure result.
    ref_ecoh = get_materials_project_B2_ecoh(chemsys)
    exp_ecoh = get_exp_B2_ecoh(chemsys)
    if exp_ecoh[0] != None:
        ax.axhline(
            y=exp_ecoh[0],
            color="royalblue",
            linestyle="-",
            label=f"{exp_ecoh[1]} (B2)",
        )
    ax.axhline(y=ref_ecoh[0], color="gray", linestyle="--", label=f"{ref_ecoh[1]} (B2)")

    ax.set_xlabel("Structure Name")
    ax.set_ylabel("Cohesive Energy (eV/atom)")
    ax.set_xticks(group_centers)
    xlabels = [x.replace("B19P", "B19'") for x in structure_names]
    ax.set_xticklabels(xlabels, rotation=0, ha="right")
    ax.tick_params(axis="x", which="both", length=0)
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim((None, -3))
    ax.legend(loc="center left", bbox_to_anchor=(1, 0.5), fontsize=6)
    # plt.legend(loc="upper right", fontsize=6)
    plt.tight_layout()
    plt.savefig(output_filename, dpi=600, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_plot = paths.figures / f"{chemsys}_CohesiveEnergyPlot.png"
    plot_cohesive_energy(paths.data / sys.argv[1], chemsys, output_plot)
