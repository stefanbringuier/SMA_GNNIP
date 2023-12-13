import numpy as np
from ase.spacegroup import crystal


def NiTi_B2_Structure(a=3.04, b=3.04, c=3.04, alpha=90.0, beta=90.0, gamma=90.0):
    """
    Austenite phase.
    """
    B2 = crystal(
        ["Ni", "Ti"],
        [(0, 0, 0), (0.5, 0.5, 0.5)],
        spacegroup=221,
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
    )

    return B2


def NiTi_B33_Structure(a=2.888, b=9.398, c=4.042, alpha=90.0, beta=90.0, gamma=90.0):
    """
    Structure from: Liu et al., Phys. Rev. B 87, 140104(R) 2018.
    Cmcm (B33/BCO) @ 0 GPa
    """
    B33 = crystal(
        ["Ni", "Ti"],
        [(0.0, 0.5847, 0.25), (0.0, 0.1440, 0.75)],
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
        spacegroup=63,
    )
    return B33


def NiTi_Cmcm_Structure(a=2.90, b=3.99, c=4.90, alpha=90.0, beta=107.01, gamma=90.0):
    """
    Wait this is B33?

    See Table IX PhysREv B 92, 134107 Lawson.
    """
    Cmcm = crystal(
        ["Ni", "Ti"],
        [(0.915362, 0.75, 0.829275), (0.641384, 0.75, 0.284766)],
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
        spacegroup=63,
    )
    return Cmcm


def NiTi_BCO_Structure(a=2.8987, b=9.3619, c=3.9870, alpha=90.0, beta=90.0, gamma=90.0):
    """
    Corresponds to materials project id: https://materialsproject.org/materials/mp-1067248
    The MP entry will show the incorrect cell, need to dowload cif to get correct cell.
    """
    BCO = crystal(
        ["Ni", "Ti"],
        [(0.0, 0.08436356, 0.75), (0.0, 0.35661795, 0.75)],
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
        spacegroup=63,
    )
    return BCO


def NiTi_B19_Structure(a=4.09, b=2.86, c=4.59, alpha=90.0, beta=90.0, gamma=90.0):
    """
    Pmma (spacegroup 51)
    https://next-gen.materialsproject.org/materials/mp-603347
    """
    B19 = crystal(
        ["Ni", "Ti"],
        [(0.25, 0.0, 0.82097441), (0.25, 0.50, 0.27440806)],
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
        spacegroup=51,
    )
    return B19


def NiTi_B19P_Structure(
    a=2.8917, b=3.9670, c=4.8256, alpha=90.0, beta=105.2296, gamma=90.0
):
    """
    I believe this is the proper B19' structure.

    From Comp. Mat Sci S. Gur paper, Seems same as https://materialsproject.org/materials/mp-1048
    """
    B19P = crystal(
        ["Ni", "Ti"],
        [(0.076622, 0.25, 0.671102), (0.369548, 0.25, 0.217074)],
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
        spacegroup=11,
    )
    return B19P


def NiTi_B19P_2_Structure(
    a=2.7998, b=4.1711, c=4.6710, alpha=90.0, beta=95.5513, gamma=90.0
):
    """
    I don't know what this structure represents, MP says its not ground state



    Structure: https://materialsproject.org/materials/mp-1179013

    The symmetrized cif file contains basis sites. Don't use the Wyckoff as
    they appear to be for the conventional cell.
    """
    B19P_2 = crystal(
        ["Ni", "Ti"],
        [(0.02095900, 0.75, 0.17804700), (0.43539200, 0.75, 0.71167800)],
        spacegroup=11,
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
    )

    return B19P_2


def NiTi_Pbcm_Structure(a=2.572, b=9.029, c=4.176, alpha=90.0, beta=90.0, gamma=90.0):
    """
    Structure from: Liu et al., Phys. Rev. B 87, 140104(R) 2018.
    Pcbm @ 20 GPa
    """

    Pbcm = crystal(
        ["Ni", "Ti"],
        [(0.5357, 0.1598, 0.25), (0.0528, 0.3975, 0.25)],
        spacegroup=57,
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
    )

    return Pbcm


def NiTi_B32_Structure(a=6.11, b=6.11, c=11.47, alpha=90.0, beta=90.0, gamma=90.0):
    """
    Structure from: Liu et al., Phys. Rev. B 87, 140104(R) 2018.
    Fd-3m (B32) @ 40 GPa

    """
    B32 = crystal(
        ["Ni", "Ti"],
        [
            (0.0, 0.0, 0.5),
            (0.5, 0.5, 0.0),
        ],
        spacegroup=227,
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
    )

    return B32


def NiTi_R_Phase_Structure(
    a=11.1606, b=11.1606, c=5.0125, alpha=90.0, beta=90.0, gamma=120.0
):
    """
    Structure: https://materialsproject.org/materials/mp-567653
    """

    R_Phase = crystal(
        ["Ni", "Ni", "Ni", "Ti"],
        [
            (0, 0, 0),
            (0.0, 0.0, 0.5),
            (0.05265937, 0.26483708, 0.20647216),
            (0.04293075, 0.21952829, 0.71750181),
        ],
        spacegroup=148,
        cellpar=[a, b, c, alpha, beta, gamma],
        pbc=True,
    )

    return R_Phase


def GetStructure(structure_name):
    structure_functions = {
        "B2": NiTi_B2_Structure,
        "B19": NiTi_B19_Structure,
        "B19P": NiTi_B19P_Structure,
        "B33": NiTi_B33_Structure,
        "Cmcm": NiTi_Cmcm_Structure,
        "BCO": NiTi_BCO_Structure,
        "Pbcm": NiTi_Pbcm_Structure,
        "B32": NiTi_B32_Structure,
        "R_Phase": NiTi_R_Phase_Structure,
    }

    if structure_name not in structure_functions:
        raise ValueError(f"Unknown structure: {structure_name}")

    structure = structure_functions[structure_name]()
    structure.info = {"structure_name": structure_name}

    return structure


def get_path_directions(path):
    """
    Calculate the principal directions along the unique segments of the band path,
    taking into account breaks indicated by commas and maintaining the sign of the directions.

    Parameters:
    path: BandPath object from ASE.

    Returns:
    List of strings representing directions in the format ["[0, 0, 0]", ...]
    """
    directions = []
    special_points = path.special_points

    # Split the path string into disconnected parts
    disconnected_paths = path.path.split(",")

    for disconnected_path in disconnected_paths:
        for i in range(len(disconnected_path) - 1):
            start_point = disconnected_path[i]
            end_point = disconnected_path[i + 1]

            start_coords = special_points[start_point]
            end_coords = special_points[end_point]

            direction = start_coords - end_coords

            direction_str = []
            for component in direction:
                if np.isclose(component, 0):
                    direction_str.append("0")
                elif component > 0:
                    direction_str.append("$\\xi$")
                else:
                    direction_str.append("$-\\xi$")
            directions.append("[{}]".format(", ".join(direction_str)))

    return directions


def spacegroup221_bandpath(
    structure,
    npoints=200,
):
    """
    The default sampling of FBZ is sufficient to use for B2.
    """
    cell = structure.cell
    path = cell.bandpath(npoints=npoints)
    special_points = path.special_points
    directions = get_path_directions(path)

    return path, directions


def spacegroup11_bandpath(
    structure,
    special_points={
        "G": [0.0, 0.0, 0.0],
        "Y": [0.5, 0.0, 0.0],
        "Z": [0.0, 0.5, 0.0],
        "B": [0.0, 0.0, 0.5],
        "C": [0.5, 0.5, 0.0],
        "A": [0.5, 0.0, 0.5],
        "D": [0.0, 0.5, 0.5],
        "E": [0.5, 0.5, 0.5],
    },
    npoints=200,
):
    """
    Ref. https://www.cryst.ehu.es/cgi-bin/cryst/programs/nph-kv-list    
    """
    cell = structure.cell
    highsym_path = "DGYAGEZG"
    path = cell.bandpath(highsym_path,special_points=special_points,npoints=npoints)
    directions = get_path_directions(path)
    return path, directions

def spacegroup51_bandpath(
    structure,
    special_points={
        "G": [0.0, 0.0, 0.0],
        "Y": [0.0, 0.5, 0.0],
        "Z": [0.0, 0.0, 0.5],
        "B": [0.0, 0.5, 0.5],
        "C": [0.5, 0.5, 0.0],
        "A": [0.5, 0.0, 0.5],
        "D": [0.5, 0.5, 0.0],
        "E": [0.5, 0.5, 0.5],
    },
    npoints=200,
):
    cell = structure.cell
    highsym_path = "BGYAGEZG"
    path = cell.bandpath(highsym_path,special_points=special_points,npoints=npoints)
    directions = get_path_directions(path)
    return path, directions


def spacegroup63_bandpath(
    structure,
    special_points={
        "G": [0.0, 0.0, 0.0],
        "Y": [-0.5, 0.5, 0.0],
        "Z": [0.0, 0.0, 0.5],
        "S": [0.0, 0.5, 0.0],
        "R": [0.0, 0.5, 0.5],
        "T": [-0.5, 0.5, 0.5],
    },
    npoints=200,
):
    """
    Taken as a < b setting
    """
    cell = structure.cell
    highsym_path = "SGZRGT"
    path = cell.bandpath(highsym_path,special_points=special_points,npoints=npoints)
    directions = get_path_directions(path)
    return path, directions


# Get a path and return the directions
def get_bandpath(structure, npoints=200):
    if structure.info["name"] == "B19P":
        path, directions = spacegroup11_bandpath(structure, npoints=npoints)
    elif structure.info["name"] == "B19":
        path, directions = spacegroup51_bandpath(structure, npoints=npoints)
    elif structure.info["name"] == "B2":
        path, directions = spacegroup221_bandpath(structure, npoints=npoints)
    elif structure.info["name"] == "BCO":
        path, directions = spacegroup63_bandpath(structure, npoints=npoints)
    else:
        raise ValueError("Not valid structure")

    return path, directions

