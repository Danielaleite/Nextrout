#!/bin/bash


local_folder=$(pwd)

base=prova_sphere
input_base=sphere_inputs.ctrl



#grid0=project_ellipsoid.dat 
#projector=project_ellipsoid.py
#folder_grid0=/home/enrico/bbmuffe/repositories/muffe_p1p0/preprocess/surface_assembly/cut_locus


folder_grid0=/home/enrico/bbmuffe/repositories/muffe_p1p0/preprocess/surface_assembly/sphere_otp
grid0=sphere_prj.dat 
projector_folder=/home/enrico/bbmuffe/repositories/muffe_p1p0/preprocess/surface_assembly/sphere_otp
projector=project_sphere.py



#sed "s|folder_grid|$folder_grid0|" inputs_sphere.ctrl > inputs.ctrl 
sed "s|folder_grid|$folder_grid0|" $input_base > inputs.ctrl 
sed -i "s|grid2use|$grid0|" inputs.ctrl


./assembly_inputs.sh $base

folder_before=${base}

for i in 1 2
do
    #
    # project previous subgrid
    # 
    cd $projector_folder
    python $projector ${local_folder}/runs/$folder_before/input/subgrid.dat ${local_folder}/runs/$folder_before/input/subgrid_prj.dat
    cd -
    #
    # set grid 
    #
    folder_grid=${local_folder}/runs/${folder_before}/input
    grid=subgrid_prj.dat
  
    sed "s|folder_grid|$folder_grid|" $input_base > inputs.ctrl 
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


