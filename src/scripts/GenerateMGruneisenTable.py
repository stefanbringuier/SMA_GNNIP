import sys
from ase.db import connect
import paths
from CalculateGruneisen import *
from ConfigsUtils import *


def GenerateMModeGruneisen(dbname, output_filename):
    db = connect(dbname)
    unique_structure_types = set(entry.structure_name for entry in db.select())

    with open(output_filename, "w") as file:
        file.write("\\begin{table}[h]\n")
        file.write("\\centering\n")

        # Begin tabular environment
        file.write("\\begin{tabular}{|c|c|}\n")
        file.write("\\hline\n")
        file.write("Structure & \\textbf{M}-Mode Gr\\\"{u}neisen Parameters \\\\\n")
        file.write("\\hline\n")

        for structure_type in unique_structure_types:
            if structure_type != "B2":
                continue
            structure_name = structure_type.replace("_", "-")
            file.write(f"{structure_name} & ")
            file.write("\\begin{tabular}{c c c c c c c}\n")  # Inner table header
            file.write("Model & TA$_1$ & TA$_2$ & LA & TO$_1$ & TO$_2$ & LO \\\\\n")
            file.write("\\hline\n")
            entries = db.select(structure_name=structure_type)
            sentries = sorted(entries, key=lambda entry: ORDER[entry.model_name])
            for entry in sentries:
                model_name =  entry.get('model_name', 'NA')
                if model_name == "ALIGNN":
                    if structure_name != "B2":
                        UserWarning(f"Skipping: {structure_name} {model_name}")
                        continue

                _, gamma = CalculateModeGruneisen(entry)
                file.write(
                    f"{model_name} & {gamma[0]:.2f} & {gamma[1]:.2f} & {gamma[2]:.2f} & {gamma[3]:.2f} & {gamma[4]:.2f} & {gamma[5]:.2f} \\\\\n"
                )

            file.write("\\hline\n")
            file.write("\\end{tabular}\\\\\n")
            file.write("\\hline\n")
        file.write("\\end{tabular}\n")
        file.write(
            "\\caption{Mode Gruneisen parameters for different NiTi structures.}\n"
        )
        file.write("\\label{tab:mode_gruneisen_niti}\n")
        file.write("\\end{table}\n")


if __name__ == "__main__":
    output_file = paths.output / "Table_NiTi_M_ModeGruneisen.tex"
    GenerateMModeGruneisen(paths.data / sys.argv[1], output_file)
