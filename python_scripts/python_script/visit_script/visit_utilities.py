#!/usr/bin/env python
import re
import os
import sys
#import argparse

def visit_utilities(sessions):
    for session in sessions:
        print session
        #if ( order == 'print'):
        command=('visit -cli -nowin -s restore_print.py '+ str(session))
        os.system(command)


if __name__ == "__main__":
    print len(sys.argv)
    if len(sys.argv) > 1:
        visit_utilities(sys.argv[1:])
    else:
        raise SystemExit("usage:  python visit_utilities.py 'order' 'session list'")

