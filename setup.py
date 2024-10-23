#!/usr/bin/env python3

import os

print('<< Installing Nextrout >>')

# Install required python packages with specific versions
python_reqs = [
    'meshpy',
    'click',
    'numpy==1.19.5',
    'scipy',
    'matplotlib==3.4.3',     
    'pandas',
    'f90wrap==0.2.7'
    'networkx',
]

for pkg in python_reqs:
    os.system('pip install ' + pkg)

# Cloning repositories
os.system('mkdir ../dmk_utilities')
repos_list = ['dmk_solver', 'globals', 'linear_algebra', 'p1galerkin', 'geometry']
for repo in repos_list:
    os.system('cd ../dmk_utilities && git clone https://gitlab.com/enrico_facca/' + repo + '.git')

## Creating builds
os.system('cd ../dmk_utilities/dmk_solver && mkdir build/ && cd build &&  cmake -DMy_Fortran_Compiler=/usr/bin/gfortran -DBLAS_DIR=/usr/lib/x86_64-linux-gnu/blas -DLAPACK_DIR=/usr/lib/x86_64-linux-gnu/lapack ../ && make')
os.system('cd ../dmk_utilities/dmk_solver/build && make f2py_interface_dmk')


print('Installation finished!')
print('In case of errors, please check the Troubleshooting section at https://gitlab.com/enrico_facca/dmk_solver')
