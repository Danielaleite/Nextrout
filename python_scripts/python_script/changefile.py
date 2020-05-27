#!/usr/bin/env python
import re
import os
import sys

def changefile(old,new,replacements):
    input = open(old)
    output = open(new,'w')

    for line in input:
        newline=line
        for rpl in replacements:
            if str(rpl[0]) in str(newline):
                temp=newline
                newline=temp.replace(str(rpl[0]), str(rpl[1]))
        output.write(newline)

    input.close()
    output.close()


if __name__ == "__main__":
    if len(sys.argv) > 1:
        changefile(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python changefile.py old new changes")
    

    
        



        
    
    
