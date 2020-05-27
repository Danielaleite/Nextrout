# README #

#
# IMPORTANT: Read READ.MD in muffe_p1p0 directory first. 
#
# This directory contains a series of scripts to 
# interface with different programs used to 
#
# - generate and pre-process inputs data
# - run simulation solving the Optimal Transport Problem
# - post-process the results
#
# Three main files controls all these programs
# 1) inputs.ctrl : (contains a series of flags and parameters to
#                   control creation of inputs data)
# 2) muffa.ctrl  : (contains a series of flags and parameters to
#                   control optimal transport solver)
# 3) dmk_folder.py ( Python script based on 
#                    "click"=Command Line Interface Creation Kit
#                    to interface with all programs)
# 
# Everything works creating and manipulating data in directory having
# -sub-folders "input"  : contains all inputs data
# -sub-folders "output" : contains all output data
# -file "muffa.fnames" : contains the list of file used by OTP solver
# 


# (STEP 0) : First use of dmk_folder script
# Use command

./dmk_folder.py --help 

# All possible commands are self explanatory.

#
# (STEP 1) : Generate inputs
# Use command

./dmk_folder assembly runs/example inputs.ctrl

# to create a directory named "example" in runs
# The flag and parameters written in inputs.ctrl 
# will control which data are generated.
# This command will call the pre-process programs
# defined in directory "muffe_p1p0/preprocess".
# Use the default  "inputs.ctrl" file as first test.
# More details in "muffe_p1p0/preprocess/READ.MD"

#
# (STEP 2) : Run Optimal Transport solver
# Use command

./dmk_folder run runs/example muffa.ctrl

# will execute the program in muffe_p1p0/code/muffe.out
# using the input data in directory "runs/example/input/".
# File "muffa.ctrl" will control the execution. It will 
# be first copied into "runs/example/input/" and then 
# the OTP solver will start working.
#
# Use the default  "muffa.ctrl" file as first test.

#
# (STEP 3) : One example to see results
# Use command

./dmk_folder optvtk runs/example muffa.ctrl

# to generate in directory "runs/example/output/vtk/"
# the files in vtk format named
# opt_tdens.vtk :: optimal transport density
# opt_pot.vtk   :: optimal Potential
# 

