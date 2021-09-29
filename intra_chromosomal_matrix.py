import cooler as cool
import numpy.lib
import os

import parameters as pars
import numpy as np


cooler = cool.Cooler("./data/GM12878-MboI-allreps-filtered.10kb.cool")
for file in os.listdir("./data/intraChromosomalMatrices"):
    if file.endswith(".txt"):
        continue
    with open(f"./data/intraChromosomalMatrices/{file}") as f, open(f"./data/intraChromosomalMatrices/{file}.txt", "w") as g:
        offset = cooler.offset(file)
        for line in f:
            line = line.rstrip("\n").split("\t")
            if len(line) != 3:
                continue
            g.write(f"{int(line[0]) - int(offset) + 1}\t{int(line[1]) - int(offset) + 1}\t{line[2]}\n")
    print(f"done with {file}")
