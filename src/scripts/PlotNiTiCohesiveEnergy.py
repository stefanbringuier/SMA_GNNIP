import sys
import matplotlib.pyplot as plt
from ase.db import connect
import paths

def GenerateCohesiveEnergyPlot(dbname, output_filename):
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
    bar_width = 0.35
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
    ax.set_title('Cohesive Energy by Structure and Model')
    ax.set_xticks([r + bar_width for r in range(len(structure_names))])
    ax.set_xticklabels(structure_names)
    ax.legend()

    # Save plot
    plt.savefig(output_filename)


if __name__ == "__main__":
    output_plot = paths.figures / "NiTi_CohesiveEnergy.png"
    GenerateCohesiveEnergyPlot(paths.data / sys.argv[1], output_plot)
