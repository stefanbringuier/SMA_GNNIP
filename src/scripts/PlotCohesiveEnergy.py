import paths
import sys
import matplotlib.pyplot as plt
from matplotlib.ticker import MaxNLocator

from ase.db import connect

from Config import ORDER

def get_materials_project_ecoh():
    """
    The obtained valued for B2 structure from  https://materialsproject.com/materials/mp-571
    
    Code is provided in this doc string as the mp_api pip package is required along with a api key.

    ```python
    with MPRester(mp_apikey) as mpr:
         docs = mpr.summary.search(material_ids=["mp-571 (B2)"])
    docs[0].energy_per_atom
    ```

    Provides the corrected energy per atom in eV.
    """
    ecoh = -7.185175
    return ecoh
    
    
def plot_cohesive_energy_test(dbname, output_filename):
    db = connect(dbname)
    unique_structure_types = set(entry.structure_name for entry in db.select())

    # Data containers
    structure_names = []
    model_names = set()
    ecoh_data = {}

    # Collecting data
    for name in unique_structure_types:
        entries = db.select(structure_name=name)
        structure_names.append(name.replace("_", "-"))
        ecoh_data[name] = {}

        for entry in entries:
            model_name = entry.get('model_name', 'NA')
            ecoh = entry.get('ecoh', 0)  # Default to 0 if not available
            ecoh_data[name][model_name] = ecoh
            model_names.add(model_name)

    structure_names = sorted(structure_names)
    model_names = sorted(model_names, key=lambda x: ORDER.get(x, float('inf')))

    # Fancy colors and patterns
    colors = ['red', 'green', 'blue', 'cyan', 'magenta', 'yellow', 'black']
    patterns = ['//', '+', '*', 'x', 'o', '||', 'O', '-', '\\']

    # Plotting
    fig, ax = plt.subplots()
    bar_width = 0.1
    space_between_bars = 0.02  # Adjust this value as needed
    indices = list(range(len(structure_names)))
    for i, model_name in enumerate(model_names):
        #offset = [(x + i * bar_width) for x in indices]
        offset = [(x + i * (bar_width + space_between_bars)) for x in indices]
        ecoh_values = [ecoh_data[struct].get(model_name, 0) for struct in structure_names]
        # Assign color and pattern
        color = colors[i % len(colors)]
        pattern = patterns[i % len(patterns)]
        ax.bar(offset, ecoh_values, bar_width, label=model_name, color=color, edgecolor='black', hatch=pattern,alpha=0.75)

    # TODO: Need to manually add ab-initio cohesive energy results.
    # Most papers are shoring difference B2-Structure in meV or
    # formation energies. For now I will use the Materials Project
    # B2 structure result.
    ref_ecoh = get_materials_project_ecoh()
    ax.axhline(y=ref_ecoh, color="gray", linestyle="--",label="mp-571 (B2)")

    ax.set_xlabel('Structure Name')
    ax.set_ylabel('Cohesive Energy (eV/atom)')
    #ax.set_title('Cohesive Energy by Structure and Model')
    ax.set_xticks([r + bar_width for r in range(len(structure_names))])
    xlabels = [ x.replace("B19P","B19'") for x in structure_names ]
    ax.set_xticklabels(xlabels, rotation=45, ha="right")
    ax.yaxis.set_major_locator(MaxNLocator(integer=True))
    ax.set_ylim((None,-3))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))

    
    plt.tight_layout()
    plt.savefig(output_filename, bbox_inches='tight')
    
def plot_cohesive_energy(dbname, output_filename):
    db = connect(dbname)
    unique_structure_types = set(entry.structure_name for entry in db.select())

    # Data containers
    structure_names = []
    model_names = set()
    ecoh_data = {}

    # Collecting data
    for structure_type in unique_structure_types:
        entries = db.select(structure_name=structure_type)
        structure_names.append(structure_type.replace("_", "-"))
        ecoh_data[structure_type] = {}

        for entry in entries:
            model_name = entry.get('model_name', 'NA')
            ecoh = entry.get('ecoh', 0)  # Default to 0 if not available
            ecoh_data[structure_type][model_name] = ecoh
            model_names.add(model_name)

    # Sorting for consistent plotting
    model_names = sorted(model_names)
    structure_names = sorted(structure_names)

    # Plotting
    fig, ax = plt.subplots()

    # Width of a bar
    bar_width = 0.1
    # Position of bars on the x-axis
    indices = list(range(len(structure_names)))

    for i, model_name in enumerate(model_names):
        # Offset each model's bar
        offset = [(x + i * bar_width) for x in indices]
        # Extract cohesive energies for each structure for the current model
        ecoh_values = [ecoh_data[struct].get(model_name, 0) for struct in structure_names]
        ax.bar(offset, ecoh_values, bar_width, label=model_name)


    # Setting the x-axis labels
    ax.set_xlabel('Structure Name')
    ax.set_ylabel('Cohesive Energy (eV/atom)')
    #ax.set_title('Cohesive Energy by Structure and Model')
    ax.set_xticks([r + bar_width for r in range(len(structure_names))])
    ax.set_xticklabels(structure_names)
    ax.legend()

    # Save plot
    plt.savefig(output_filename)


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_plot = paths.figures / f"{chemsys}_CohesiveEnergyPlot.png"
    plot_cohesive_energy_test(paths.data / sys.argv[1], output_plot)
