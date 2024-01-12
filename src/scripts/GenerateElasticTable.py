import sys
from ase.db import connect
import paths

from TableConfig import SPACEGROUP_MAP, ORDER


def generate_elastic_constants_table(dbname, chemsys, output_filename):
    db = connect(dbname)
    
    sorder = ORDER[chemsys]["structure"]
    default_order = max(sorder.values())
    unique_structure_types = set(entry.structure_name for entry in db.select())
    structure_order = sorted(unique_structure_types, key=lambda x: sorder.get(x, default_order))


    with open(output_filename, "w") as file:
        # Top-level table with two main columns
        file.write("\\begin{longtable}{|l|l|}\n")
        file.write(
            "\\caption{Elastic constants for "
            + chemsys
            + ".} \\label{tab:elastic_"
            + chemsys
            + "}\\\\\n"
        )
        file.write("\\hline\n")
        file.write("Structure (\#) & Properties \\\\\n")
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

        for structure_type in structure_order:
            structure_name = structure_type.replace("_", "-")
            spg_num = SPACEGROUP_MAP[structure_name]
            file.write(f"{structure_name} ({spg_num}) & ")

            entries = db.select(chemsys=chemsys, structure_name=structure_type)

            # Filter and sort based on models defined in ORDER
            sentries = sorted(
                [
                    entry
                    for entry in entries
                    if entry.model_name in ORDER[chemsys]["model"]
                ],
                key=lambda entry: ORDER[chemsys]["model"][entry.model_name],
            )

            file.write("\\begin{tabularx}{\columnwidth}{X X X X X X X X X X X X X X }\n")
            file.write(
                "Model & C$_{11}$ & C$_{22}$ & C$_{33}$ & C$_{12}$ & C$_{13}$ & C$_{23}$ & C$_{44}$ & C$_{55}$ & C$_{66}$ & C$_{16}$ & C$_{26}$ & C$_{36}$ & C$_{45}$ \\\\\n"
            )
            file.write("\\hline\n")

            for entry in sentries:
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
                    val = elastic_data.get(ec, "-")
                    sval = f"{val:0.0f}" if val != "-" else val
                    if i < len(c) - 1:
                        file.write(f"{sval} & ")
                    else:
                        file.write(f"{sval} \\\\\n")
                        
            file.write("\\end{tabularx} \\\\\n")
            file.write("\\hline\n")

        file.write("\\end{longtable}\n")

    return None

if __name__ == "__main__":
    dbname = sys.argv[1]
    chemsys = sys.argv[2]
    output_file = paths.output / f"Table_{chemsys}_Elastic_Constants.tex"
    generate_elastic_constants_table(paths.data / dbname, chemsys, output_file)
