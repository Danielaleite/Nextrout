#!/bin/bash

base=conv_sphere
folder_grid0=/home/enrico/bbmuffe/repositories/muffe_p1p0/preprocess/surface_assembly/sphere_otp
grid0=sphere_prj.dat 

sed "s|folder_grid|$folder_grid0|" inputs_sphere.ctrl > inputs.ctrl 
sed -i "s|grid2use|$grid0|" inputs.ctrl


./assembly_inputs.sh $base

folder_before=${base}

for i in 1 2 3 4
do
    #
    # project previous subgrid
    # 
    echo 'before' $folder_before
    cd ${folder_grid0}
    echo 'cur dir'
    pwd
    python project_sphere.py /home/enrico/bbmuffe/dmkp1p0/runs/$folder_before/input/subgrid.dat /home/enrico/bbmuffe/dmkp1p0/runs/$folder_before/input/subgrid_prj.dat
    cd -
      
    #
    # set grid 
    #
    folder_grid=/home/enrico/bbmuffe/dmkp1p0/runs/${folder_before}/input
    grid=subgrid_prj.dat
  
    sed "s|folder_grid|$folder_grid|" inputs_sphere.ctrl > inputs.ctrl 
    sed -i "s|grid2use|$grid|" inputs.ctrl
    
    #
    # set folder name 
    #
    folder=${base}_nref${i}
    echo $folder

    #
    # create inputs
    #
    ./assembly_inputs.sh $folder

    #
    # set name for next start
    #
    folder_before=$folder
done


