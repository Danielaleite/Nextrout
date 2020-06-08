#!/usr/bin/env python
import numpy as np
import sys
import os
import timedata as td

fsource=sys.argv[1]
fsink=sys.argv[2]
fforcing=sys.argv[3]
try:
   option=sys.argv[4]
except:
   option='source'


source=td.read_steady_timedata(fsource)
sink=td.read_steady_timedata(fsink)

if ( len(source) !=  len(sink) ):
   print ('Dimension mismatch sink and source')

forcing=np.zeros(len(source))

if (option == 'source'):
   for i in range(len(source)):
      forcing[i]=source[i]
      if ( forcing[i] == 0.0 ) :
         forcing[i]=-sink[i]

if (option == 'sink'):
   for i in range(len(source)):
      forcing[i]=-sink[i]
      if ( forcing[i] == 0.0 ) :
         forcing[i]=source[i]

if (option == 'diff'):
   for i in range(len(source)):
      forcing[i]=source[i]-sink[i]

td.write_steady_timedata(fforcing,forcing)
