#!usr/bin/env python3

import os

#TODO: add python requirements as well

# following enricos instructions

## first move to the dmk_utilities

## then clone every

repos_list = []#['dmk_solver','globals', 'linear_algebra','p1galerkin','geometry']
for repo in repos_list: os.system('cd dmk_utilities && git clone https://gitlab.com/enrico_facca/'+repo+'.git')

## create build/'s


os.system('cd dmk_utilities/dmk_solver && mkdir build/ &&  cd build &&  cmake ../ &&  make')
os.system('cd dmk_utilities/dmk_solver/build && make f2py_interface_dmk')


# then care abt this part

'''

move_folders = {
    "otp_utilities/par_files": [
        "cp dmk_folder_structure.py ../../python_scripts/python_script"
    ],
    "otp_utilities/par_files": ["cp utilities.py ../../python_scripts/python_script"],
}
execute(move_folders)

move_folders2 = {
    "otp_utilities/": [" cp -r par_files ./muffe_sparse_optimization/simplifications/"],
}
execute(move_folders2)

move_folders3 = {
    "otp_utilities/par_files": [
        " cp muffa.fnames ../muffe_sparse_optimization/simplifications/"
    ],
}
execute(move_folders3)

move_folders4 = {
    "otp_utilities/par_files": [
        " cp  dmk_folder_structure.py ../muffe_sparse_optimization/simplifications/python_script"
    ],
}
execute(move_folders4)


move_folders5 = {
    "otp_utilities/par_files": [
        "cp utilities.py ../muffe_sparse_optimization/simplifications/python_script"
    ],
}
execute(move_folders5)

make_dir = {
    ".": [" mkdir results", "mkdir results/discrete", "mkdir results/continuous"]
}
execute(make_dir)
'''