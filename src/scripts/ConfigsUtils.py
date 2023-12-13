ORDER = {
    "Mutter": 1,
    "Zhong": 2,
    "M3GNet": 3,
    "CHGNet": 4,
    "MACE": 5,
    "ALIGNN": 6,
}


SPACEGROUP_MAP = {"B2" : 221,
                  "B19P" : 11,
                  "B19" : 11,
                  "BCO" : 63,
                  }

def GetPhononCalcConfig(structure_name, potential):
    """
    Settings have been manually calibrated.
    """
    # Supercell, displacement (ang.), band path points
    configurations = {
        "B2": {
            "Mutter": ((8, 8, 8), 0.0935, 200),
            "Zhong": ((8, 8, 8), 0.0935, 200),
            "M3GNet": ((8, 8, 8), 0.260, 200),
            "CHGNet": ((8, 8, 8), 0.030, 200),
            "MACE": ((8, 8, 8), 0.010, 200),
            "ALIGNN": ((6, 6, 6), 0.050, 200),
        },
        "B19": {
            "Mutter": ((7, 9, 7), 0.05, 200),
            "Zhong": ((7, 9, 7), 0.05, 200),
            "M3GNet": ((7, 9, 7), 0.05, 200),
            "CHGNet": ((7, 9, 7), 0.05, 200),
            "MACE": ((7, 9, 7), 0.05, 200),
            "ALIGNN": ((5, 7, 5), 0.05, 200),
        },
        "B19P": {
            "Mutter": ((8, 6, 4), 0.05, 200),
            "Zhong": ((8, 6, 4), 0.05, 200),
            "M3GNet": ((8, 6, 4), 0.05, 200),
            "CHGNet": ((8, 6, 4), 0.05, 200),
            "MACE": ((8, 6, 4), 0.05, 200),
            "ALIGNN": ((7, 5, 3), 0.05, 200),
        },
        "BCO": {
            "Mutter": ((8, 4, 6), 0.05, 200),
            "Zhong": ((8, 4, 6), 0.05, 200),
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
