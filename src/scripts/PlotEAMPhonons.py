import sys
import paths
import json
from PlotPhonons import *


eam_models = ["Mutter", "Zhong"]
dbname = paths.data / sys.argv[1]
nstrains = int(sys.argv[2])
structures = sys.argv[3:]

for m in eam_models:
    for s in structures:
        PlotPhonons(dbname, m, s, num_strains=nstrains)
