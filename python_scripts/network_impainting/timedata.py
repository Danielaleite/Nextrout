#!/usr/bin/env python
import numpy as np
import sys
import os


def read_steady_timedata(filename):
   fin=open(str(filename), 'r')
   lines = fin.readlines()
   dimdata=1
   ndata=int(lines[0].split()[1])
   data=np.zeros(ndata)   
   ninputs= int(lines[2])
   for i in range(ninputs):
      line=lines[3+i]
      idata=int(line.split()[0])-1
      data[idata] = float(line.split()[1])
   fin.close()
   return data;

def write_steady_timedata(filename,data):
   fout=open(str(filename), 'w')
   dimdata=1
   ndata=len(data)
   ninputs = np.count_nonzero(data)
   print ninputs
   fout.write(str(dimdata)+" "+ str(ndata)+"\n")
   fout.write("time  0.0 \n")
   fout.write(str(ninputs)+" \n")
   for i in range(ndata):
      if (data[i] !=0.0 ):
         fout.write(str(i+1) + ' ' +str(data[i])+"\n")
   fout.write("time  1.0e30 \n")
   fout.close()
   return;

def balance(data,option=None):
    if (option==None):
        option='plus'
    plus=0.0
    minus=0.0
    
    balanced_data=np.zeros(len(data))

    for i in range(len(data)):
        if (data[i] > 0.0) :
            plus=plus+data[i]
        else:
            minus=minus-data[i]

    if (option == 'plus'):
        print (str(option)+'Positive part preserved | Negative part scaled')
        for i in range(len(data)):
            if ( data[i] > 0.0 ):
               balanced_data[i] = data[i]
            if ( data[i] < 0.0 ) :
               balanced_data[i] = data[i]*plus/minus

    if (option == 'minus'):
        print (str(option)+'Positive part scaled | Negative part preserved')
        for i in range(len(data)):
            if ( data[i] > 0.0 ):
                balanced_data[i] = data[i]*minus/plus
            if ( data[i] < 0.0 ) :
                balanced_data[i] = data[i]

    print ( 'Mass positive part:=',plus)
    print ( 'Mass negative part:=',minus)
    print ( 'Input  unbalanced :=',data.sum())
    print ( 'Output unbalanced:=',balanced_data.sum())

    return balanced_data;

