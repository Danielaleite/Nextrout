#!/bin/bash
#
# program to convert all corrupted netwotk images to file  .dat and .vtk
#


folder_rep=$(head ../location_repository.txt)

example_folder=$1
example_folder=$(readlink -f $example_folder)
factor=$2


#
# check input existence of other data
#
for folder in corropted_network
do
    if [ ! -d $example_folder/$folder ] 
    then
	echo ${folder}' does not exist'
    fi 
done

folder_data=${example_folder}/common_data/
folder_corr=${example_folder}/corrupted_data/
mkdir ${folder_corr}

#
# convert all corropted data into 
#



for img in ${example_folder}/corrupted_network/*.png
do
    echo $img
    fname=$(basename ${img})
    clean_name=$(echo ${fname%.*})
    data_name=${clean_name}'.dat'
    echo $clean_name
    python image2dat.py ${img} ${folder_corr}/$data_name ${factor} 
    # copy required because fortran can not open in different logical unit the same files
    cp ${folder_corr}/$data_name ${folder_corr}/bis_$data_name
    #
    # convert to vtk
    #
    $folder_rep/../geometry/timedata2vtk/timedata2vtk.out ${folder_data}/grid.dat ${folder_corr}/$data_name ${folder_corr}/
    mv ${folder_corr}/${clean_name}00*.vtk ${folder_corr}/${clean_name}.vtk 
done
