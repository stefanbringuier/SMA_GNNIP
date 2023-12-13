import paths
from ase.db import connect


def set_json_db_metadata(dbname):
    """
    Sets metadata for a JSON style ASE database.

    Arguments:
        dbname (string): full path to JSON ASE database.

    Returns:
        None

    """
    db = connect(dbname, append=True)
    db.metadata = {
        "title": "NiTi Calculations",
        "key_descriptions": {
            "vi": ("Initial Volume", "Initial volume of structure", "Å^3"),
            "v0": ("Equilibrium Volume", "Equilibrium volume from EOS fit", "Å^3"),
            "e0": ("Equilibrium Energy", "Equilibrium energy from EOS fit", "eV"),
            "bulk_modulus": ("Bulk Modulus", "Bulk modulus from EOS fit", "GPa"),
            "fit_eq": ("Fit Equation", "Equation used for the EOS fit", "string"),
            "strain": ("Strain", "Isotropic Strain", "fractional value"),
            "potential": ("Potential", "Potential used in the calculation", "string"),
            "spacegroup": ("Spacgroup #", "International spacegroup number", "number"),
            "ecoh": ("Cohesive energy", "Cohesive energy per atom", "eV/atom"),
            "structure_name": ("Structure Name", "Prototype structure name", "string"),
            "model_name": ("Model Name", "Name of the model/parameters", "string"),
        },
        "default_columns": [
            "id",
            "calculator",
            "formula",
            "potential",
            "energy",
            "vi",
            "v0",
            "e0",
            "ecoh",
            "bulk_modulus",
            "fit_eq",
            "strain",
            "volume",
            "structure_name",
            "model_name",
            "spacegroup",
        ],
    }
    return None


if __name__ == "__main__":
    db_path_file = paths.data / "NiTi_Structures.json"
    set_json_db_metadata(db_path_file)
