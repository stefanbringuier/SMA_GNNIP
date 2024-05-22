import matplotlib.pyplot as plt
import paths
from Structures import *

structure_names = ["B2", "B19", "B19P", "BCO"]
fmap = {
    "B2": NiTi_B2_Structure,
    "B19": NiTi_B19_Structure,
    "B19P": NiTi_B19P_Structure,
    "BCO": NiTi_BCO_Structure,
}

for name in structure_names:
    structure = fmap[name]()
    structure.info["name"] = name
    spg = structure.info["spacegroup"].no
    path, _ = get_bandpath(structure, npoints=1)
    bz_diagram = f"{name}_BrillouinZonePointsSampled.png"
    bz_plot = path.plot(vectors=True)
    bz_plot.figure.savefig(paths.figures / bz_diagram)
    special_points = path.special_points

    # Write table of special points
    latex_table = "\\begin{table}[h]\n"
    latex_table += "\\centering\n"
    latex_table += "\\begin{tabular}{|c c|}\n"
    latex_table += "\\hline\n"
    latex_table += "Point & Coordinates \\\\"
    latex_table += "\\hline\n"

    for point, coords in special_points.items():
        coords_str = ", ".join(map(str, coords))
        latex_table += f"{point} & {coords_str} \\\\\n"

    latex_table += "\\hline\n"
    latex_table += "\\end{tabular}\n"
    latex_table += f"\\caption{{Special points sampled in the irreducible Brillouin Zone for space group {spg}}}\n"
    latex_table += "\\end{table}"

    sp_table = f"{name}_SpecialSymmetryPointsBZ.tex"
    with open(paths.output / sp_table, "w") as file:
        file.write(latex_table)
