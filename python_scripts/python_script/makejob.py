7#!/usr/bin/env python
import sys
import os
from changefile import changefile 

def makejob(base,launch_folder,jobname,folder,longshort):
    changes=[];
    changes.append(['path_folder',str(folder)])
    changes.append(['jobname',str(jobname)])
    changes.append(['local_folder',str(launch_folder)])
    if ( longshort == 'long' ) :
        changes.append(['cluster_type','cluster_long'])
        changes.append(['maxtime','999:00:00'])

    elif ( longshort == 'short' ) :
        changes.append(['cluster_type','cluster_short'])
        changes.append(['maxtime','00:59:00'])

    
    newjobname=str(launch_folder)+'/jobs/'+jobname+'.job'
    
    changefile(base,newjobname,changes)

if __name__ == "__main__":
    if len(sys.argv) > 1:
        makejob(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python makejob.py base new changes")
