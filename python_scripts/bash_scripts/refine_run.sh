#!/bin/bash

#########################################
# definitions
repo_path=$(head location_repository.txt)
echo $repo_path
muffa=${repo_path}/code/muffa.out
controls_file=muffa.ctrl
REFINE_DATA=${repo_path}/../geometry/refine_data/code/refine_data.out

fnames=muffa.fnames
controls_file=muffa.ctrl
###########################################

###########################################
# read inputs folder and flags
folder=$1
nref=$2
var_ref=$3

echo $folder
echo $nref
echo $var_ref



# copy controls files
cp $controls_file $folder/input
sed -i "32s|.*|$var_ref                ! var_tdens|" $folder/input/muffa.ctrl

i=0
cd $folder
pwd
echo "local"
$muffa
cd -

for i in `seq 1 $nref`
do 
    # folders name
    folder_old=$folder
    folder=${folder}_nref${i}

    #get last time
    last_time=$(tail -1 ${folder_old}/output/timefun/time.dat)
     
    # use final tdens as new initial data
    $REFINE_DATA\
    ${folder_old}/output/result/opt_tdens.dat \
    ${folder_old}/input/subgrid.dat \
    ${folder_old}/input/parent.dat \
    ${folder}/input/tdens0.dat
   
       
    cp muffa.ctrl $folder/input/ 
    # shift time set not time evolution 
    sed -i "19s|.*|$last_time                  ! tzero|" ${folder}/input/muffa.ctrl
    if [ $i != $nref ];
    then
	sed -i "32s|.*|$var_ref                  ! var_tdens|" ${folder}/input/muffa.ctrl
    fi
   
    # run program
    cd $folder
    $muffa
    cd -    
done

