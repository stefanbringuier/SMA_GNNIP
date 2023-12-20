ORDER = {
    "Mutter": 1,
    "Zhong": 2,
    "Ko" : 3,
    "M3GNet": 4,
    "CHGNet": 5,
    "MACE": 6,
}

SUBPLOT_ORDER = {
    "a": "Mutter",
    "b": "Zhong",
    "c": "Ko",
    "d": "M3GNet",
    "e": "CHGNet",
    "f": "MACE",
}


SPACEGROUP_MAP = {"B2" : 221,
                  "B19P" : 11,
                  "B19" : 51,
                  "BCO" : 63,
                  }

def get_phonon_config(structure_name, potential):
    """Settings for the phonon calculations

    - supercell
    - atomic displacement (angstrom)
    - q-points along band path

    Settings have been manually calibrated.
    """
    # Supercell, displacement (ang.), band path points
    configurations = {
        "B2": {
            "Mutter": ((8, 8, 8), 0.0935, 200),
            "Zhong": ((8, 8, 8), 0.0935, 200),
            "Ko": ((8, 8, 8), 0.01, 200),
            "M3GNet": ((8, 8, 8), 0.260, 200),
            "CHGNet": ((8, 8, 8), 0.030, 200),
            "MACE": ((8, 8, 8), 0.010, 200),
            "ALIGNN": ((6, 6, 6), 0.050, 200),
        },
        "B19": {
            "Mutter": ((7, 9, 7), 0.05, 200),
            "Zhong": ((7, 9, 7), 0.05, 200),
            "Ko": ((7, 9, 7), 0.03, 200),
            "M3GNet": ((7, 9, 7), 0.05, 200),
            "CHGNet": ((7, 9, 7), 0.05, 200),
            "MACE": ((7, 9, 7), 0.05, 200),
            "ALIGNN": ((5, 7, 5), 0.05, 200),
        },
        "B19P": {
            "Mutter": ((8, 6, 4), 0.05, 200),
            "Zhong": ((8, 6, 4), 0.05, 200),
            "Ko": ((8, 6, 4), 0.01, 200),
            "M3GNet": ((8, 6, 4), 0.05, 200),
            "CHGNet": ((8, 6, 4), 0.05, 200),
            "MACE": ((8, 6, 4), 0.05, 200),
            "ALIGNN": ((7, 5, 3), 0.05, 200),
        },
        "BCO": {
            "Mutter": ((8, 4, 6), 0.05, 200),
            "Zhong": ((8, 4, 6), 0.05, 200),
            "Ko": ((8, 4, 6), 0.05, 200),
            "M3GNet": ((8, 4, 6), 0.05, 200),
            "CHGNet": ((8, 4, 6), 0.05, 200),
            "MACE": ((8, 4, 6), 0.05, 200),
            "ALIGNN": ((5, 2, 3), 0.05, 200),
        },
    }

    # Validate structure_name and potential_key
    if structure_name not in configurations or potential not in configurations[structure_name]:
        raise ValueError(f"Configuration not found for structure '{structure_name}' and potential '{potential_key}'")

    # Retrieve configuration
    return configurations[structure_name][potential]
