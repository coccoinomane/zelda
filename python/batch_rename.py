#!/usr/bin/env python
import os
import sys
import re

folder = sys.argv[1]

n_renamed = 0

for file in os.listdir(folder):
    # if "z989" in file:
    #     os.rename(os.path.join(folder,file), re.sub("z989", "z0.989", os.path.join(folder,file)))
    #     n_renamed = n_renamed+1
    if "ngrid1." in file:
        if "isolation" in file:
            os.rename(os.path.join(folder,file), re.sub("ngrid1.", "ngrid12_neigh1.", os.path.join(folder,file)))
            n_renamed = n_renamed+1
        
print "Renamed %d files" % n_renamed