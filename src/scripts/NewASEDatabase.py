import sys
import paths
from ase.db import connect
from SetMetadata import set_db_metadata

dbname = sys.argv[1]
db = connect(paths.data / dbname, append=False)

#Need to write something otherwise it doesn't work
from ase import Atoms
from ase.calculators.emt import EMT
h2 = Atoms('H2', [(0, 0, 0), (0, 0, 0.7)])
h2.calc = EMT()
h2.get_total_energy()
h2.get_forces()
db.write(h2)
db.delete(ids=[1])
set_db_metadata(paths.data / dbname)

exit()
