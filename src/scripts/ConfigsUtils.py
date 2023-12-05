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
            "ALIGNN": ((6, 6, 6), 0.025, 200),
            "Mutter": ((4, 4, 4), 0.025, 200),
            "Zhong": ((5, 5, 5), 0.025, 200),
            "M3GNet": ((5, 5, 5), 0.01, 200),
            "CHGNet": ((5, 5, 5), 0.050, 200),
            "MACE": ((5, 5, 5), 0.01, 200),
        },
        "B19": {
            "Mutter": ((4, 4, 4), 0.050, 200),
            "Zhong": ((5, 5, 5), 0.050, 200),
            "M3GNet": ((4, 4, 4), 0.010, 200),
            "CHGNet": ((5, 5, 5), 0.010, 200),
            "MACE": ((4, 4, 4), 0.01, 200),
            "ALIGNN": ((4, 4, 4), 0.025, 200),
        },
        "B19P": {
            "Mutter": ((4, 4, 4), 0.050, 200),
            "Zhong": ((5, 5, 5), 0.050, 200),
            "M3GNet": ((4, 4, 4), 0.010, 200),
            "CHGNet": ((5, 5, 5), 0.010, 200),
            "MACE": ((4, 4, 4), 0.01, 200),
            "ALIGNN": ((4, 4, 4), 0.025, 200),
        },
        "BCO": {
            "Mutter": ((2, 2, 2), 0.050, 200),
            "Zhong": ((2, 2, 2), 0.050, 200),
            "M3GNet": ((2, 2, 2), 0.010, 200),
            "CHGNet": ((2, 2, 2), 0.010, 200),
            "MACE": ((2, 2, 2), 0.01, 200),
            "ALIGNN": ((2, 2, 2), 0.025, 200),
        },
    }

    # Validate structure_name and potential_key
    if structure_name not in configurations or potential not in configurations[structure_name]:
        raise ValueError(f"Configuration not found for structure '{structure_name}' and potential '{potential_key}'")

    # Retrieve configuration
    return configurations[structure_name][potential]
