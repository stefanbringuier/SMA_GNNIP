import sys
from ase.db import connect
import paths


def generate_elastic_constants_table(dbname, chemsys, output_filename):
    db = connect(dbname)

    with open(output_filename, "w") as file:
        # Top-level table with two main columns
        file.write("\\begin{longtable}{|l|c|}\n")
        file.write(
            "\\caption{Elastic constants for "
            + chemsys
            + ".} \\label{tab:elastic_"
            + chemsys
            + "}\\\\\n"
        )
        file.write("\\hline\n")
        file.write("Structure & Properties \\\\\n")
        file.write("\\hline\n")
        file.write("\\endfirsthead\n")
        file.write(
            "\\caption[]{Elastic constants for " + chemsys + ". (Continued)}\\\\\n"
        )
        file.write("\\hline\n")
        file.write("Structure & Properties \\\\\n")
        file.write("\\endhead\n")
        file.write("\\hline\n")
        file.write("\\endfoot\n")

        # Group entries by structure type
        grouped_entries = {}
        for entry in db.select(chemsys=chemsys):
            structure_name = entry.structure_name.replace("_", "-")
            if structure_name not in grouped_entries:
                grouped_entries[structure_name] = []
            grouped_entries[structure_name].append(entry)

        # Process each group
        for structure_name, entries in grouped_entries.items():
            file.write(f"{structure_name} & ")

            # Subtable for model and elastic constants
            file.write("\\begin{tabular}{c|c|c|c|c|c|c|c|c|c|c|c|c|c}\n")
            file.write(
                "Model & C$_{11}$ & C$_{22}$ & C$_{33}$ & C$_{12}$ & C$_{13}$ & C$_{23}$ & C$_{44}$ & C$_{55}$ & C$_{66}$ & C$_{16}$ & C$_{26}$ & C$_{36}$ & C$_{45}$ \\\\\n"
            )
            file.write("\\hline\n")

            for entry in entries:
                model_name = entry.get("model_name", "NA")
                elastic_data = entry.data.get("elastic_constants", {})
                file.write(f"{model_name} & ")
                c = [
                    "C_11",
                    "C_22",
                    "C_33",
                    "C_12",
                    "C_13",
                    "C_23",
                    "C_44",
                    "C_55",
                    "C_66",
                    "C_16",
                    "C_26",
                    "C_36",
                    "C_45",
                ]
                for i, ec in enumerate(c):
                    value = elastic_data.get(ec, "-")
                    if i < len(c):
                        file.write(f"{value} & ")
                file.write("\\\\\n")

            file.write("\\end{tabular} \\\\\n")
            file.write("\\hline\n")

        file.write("\\end{longtable}\n")

    return None


if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_file = paths.output / f"Table_{chemsys}_Elastic_Constants.tex"
    generate_elastic_constants_table(paths.data / dbname, chemsys, output_file)
