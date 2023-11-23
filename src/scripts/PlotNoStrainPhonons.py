import sys
import paths
from PlotPhonons import *


models = ["Mutter", "Zhong","M3GNet","CHGNet","MACE","ALIGNN"]
dbname = paths.data / sys.argv[1]
structures = sys.argv[2:]

for s in structures:
    PlotAllModelPhonons(dbname, models, s, strain=0.00)
