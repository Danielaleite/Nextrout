#!/usr/bin/env python
from PIL import Image
import numpy as np
import sys
import os
import subprocess


def image2dat(img_name,ndisk,radius,mask_name=None,factor=100):
   ndisk=int(ndisk)
   radius=int(radius)
   if ( factor != '100' ):
      cmd='convert -resize '+str(factor)+'%  '+ str(img_name)+'  img_temp'
      subprocess.call(cmd,shell=True)
      img_name='img_temp'

   #open file in fileList:
   img_file = Image.open(img_name)


   # get original image parameters...
   height, width = img_file.size

   
   data_mask=np.zeros([width,height])
   data_mask[:,:]=255

   check = 2
   
   for i in range(ndisk) :
      x_center=int(np.random.rand()*width)
      y_center=int(np.random.rand()*height)
      
      for ip in range( x_center - radius - check,  x_center + radius + check):
         for jp in range( y_center - radius - check,  y_center + radius + check) :
            if ( ip <width ) and ( jp < height ) :
               if ( (ip-x_center)**2 +(jp-y_center)**2 < radius**2 ):
                  data_mask[ip,jp]=0
               
   print 
   img = Image.fromarray(data_mask)
   img = img.convert("L")
   img.save('mask.png')


if __name__ == "__main__":
    if len(sys.argv) > 1:
       image2dat(*sys.argv[1:])
    else:
        raise SystemExit("usage:  python  image2dat image data [grid]")
