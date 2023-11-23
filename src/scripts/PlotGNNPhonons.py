import sys
import paths
import json
from PlotPhonons import *


gnn_models = ["M3GNet","CHGNet","MACE"]
dbname = paths.data / sys.argv[1]
nstrains = int(sys.argv[2])
structures = sys.argv[3:]

for m in gnn_models:
    for s in structures:
        PlotPhonons(dbname, m, s, num_strains=nstrains)
