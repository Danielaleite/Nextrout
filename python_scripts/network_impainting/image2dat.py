#!/usr/bin/env python
from PIL import Image
import numpy as np
import sys
import os
import subprocess

import rectangle_grid as grid 

def createFileList(myDir, format='.jpg'):
   fileList = []
   print(myDir)
   for root, dirs, files in os.walk(myDir, topdown=False):
      for name in files:
         if name.endswith(format):
            fullName = os.path.join(root, name)
            fileList.append(fullName)
   return fileList

# load the original image
#myFileList = createFileList('path/to/directory/')

def image2dat(img_name,dat_name,factor=100,grid_name=None):
   if ( factor != '100' ):
      cmd='convert -resize '+str(factor)+'%  '+ str(img_name)+'  img_temp'
      subprocess.call(cmd,shell=True)
      img_name='img_temp'

   #open file in fileList:
   img_file = Image.open(img_name)
   # img_file.show()

   # get original image parameters...
   width, height = img_file.size
   #compress
   # if ( factor !=1.0):
   #    bsize = int(float(img_file.size[0]) * float(factor))
   #    hsize = int(float(img_file.size[1]) * float(factor))
   #    img_file = img_file.resize((bsize, hsize), Image.ANTIALIAS)


   format = img_file.format
   mode = img_file.mode

   # Make image Greyscale
   img_grey = img_file.convert('L')
   #img_grey.save('result.png')
   #img_grey.show()

   # Save Greyscale values
   print ('Pixels : '+str(img_grey.size[1])+' x '+ str(img_grey.size[0]))
   value = np.asarray(img_grey.getdata(), dtype=np.float)
   value[:]=value[:]/255.0
   value[value>0]=1

   
   value=value.reshape((img_grey.size[0], img_grey.size[1]))
   #value = value.flatten()
   img_grey.save('result.png')
   #print(value)
   ntdens = 2*value.shape[0]*value.shape[1]
   ninputs = 2*np.count_nonzero(value)

   
   f_out=open(str(dat_name), 'w')
   f_out.write("1 "+ str(ntdens)+" \n")
   f_out.write("time  0.0 \n")
   f_out.write(str(ninputs)+" \n")
   n=0
   for i in range(value.shape[0]):
      for j in range(value.shape[1]):
         n=n+1
         if (value[i,j] != 0.0):
            f_out.write(str(n) + ' ' +str(float(value[i,j]))+"\n")
         #f_out.write(str(float(value[i,j]))+"\n")
         n=n+1
         if (value[i,j] != 0.0):
            f_out.write(str(n) + ' ' +str(float(value[i,j]))+"\n")
   f_out.write("time  1.0e30 \n")
   f_out.close()


   if ( grid_name!=None):
      grid.rectangule_grid(value.shape[0],value.shape[1],grid_name)   
   return;


if __name__ == "__main__":
    if len(sys.argv) > 1:
       image2dat(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python  image2dat image data [grid]")
