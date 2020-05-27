#!/bin/bash
#
# scripts to prepare all common inputs used for different penalties
# in creates inputs file grid.dat subgrid.dat rhs_integrated.dat
#
example_folder=$1
echo $example_folder
example_folder=$(readlink -f ${example_folder})
echo $example_folder
factor=$2
if [ -z $factor ]
then
    factor=100
fi

#
# check input images
#
for file in source.png sink.png network.png
do
    if [ ! -f $example_folder/$file ] 
    then
	echo ${file}' not found'
    fi 
done



#
# prepapre source sink data and grid
#
folder_data=$example_folder/common_data
mkdir ${folder_data}

#
# convert images to data
#
python image2dat.py $example_folder/source.png ${folder_data}/source.dat ${factor} ${folder_data}/grid.dat
python image2dat.py $example_folder/sink.png ${folder_data}/sink.dat   ${factor} 



#
# create balances source and sink term
#
python source_plus_sink_balance.py ${folder_data}/source.dat ${folder_data}/sink.dat $folder_data/forcing.dat
python forcing2sourcesink.py ${folder_data}/forcing.dat ${folder_data}/source.dat $folder_data/sink.dat
rm $folder_data/forcing.dat


#
# convert 2 data real network
# 
python image2dat.py ${example_folder}/network.png ${folder_data}/network.dat ${factor} 

#
#
#
folder_rep=$(head ../location_repository.txt)
$folder_rep/../geometry/timedata2vtk/timedata2vtk.out ${folder_data}/grid.dat ${folder_data}/network.dat ${folder_data}/
mv ${folder_data}/network00*.vtk ${folder_data}/network.vtk 


#
# prepar otp common data
#
abspath=$(readlink -f ${folder_data})

sed "s|keyword|$abspath|g" default.ctrl>../inputs.ctrl

cd ..
./dmk_folder.py assembly ${folder_data}/opt inputs.ctrl
cd -
