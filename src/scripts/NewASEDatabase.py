import sys

import paths
from ase import Atoms
from ase.calculators.emt import EMT
from ase.db import connect


def new_ase_database(dbname):
    """
    Create a ASE database with a dummy entry then delete.

    Arguments:
        dbname (string): full path to ASE database file to be created.

    Returns:
        None

    Notes:
        A dummy entry is created then deleted so the index will
        be shifted by 1.

    """
    db = connect(paths.data / dbname, append=False)

    # NOTE: Must write something to the connected database
    # but then delete it.
    h2 = Atoms("H2", [(0, 0, 0), (0, 0, 0.7)])
    h2.calc = EMT()
    h2.get_total_energy()
    h2.get_forces()
    db.write(h2)
    db.delete(ids=[1])

    return None


if __name__ == "__main__":
    new_ase_database(sys.argv[1])
