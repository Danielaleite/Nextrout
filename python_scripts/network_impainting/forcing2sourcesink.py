#!pyhton
import numpy as np
import sys
import timedata as td

fforcing=sys.argv[1]
fsource=sys.argv[2]
fsink=sys.argv[3]

forcing=td.read_steady_timedata(fforcing)

source=np.zeros(len(forcing))
sink=np.zeros(len(forcing))

for i in range(len(forcing)):
    if ( forcing[i]>0.0 ) :
        source[i]=forcing[i]
    else:
        sink[i]=-forcing[i]

td.write_steady_timedata(fsource,source)
td.write_steady_timedata(fsink,sink)
