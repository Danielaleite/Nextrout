#!/usr/bin/env python3

import os

print('<< Installing Nextrout >>')

# Create a virtual environment
os.system('python3 -m venv nextrout-env')
os.system('source nextrout-env/bin/activate')

# Upgrade pip to the latest version
os.system('pip install --upgrade pip')

# Install required python packages with specific versions
python_reqs = [
    'meshpy',
    'click',
    'numpy==1.22.0',          
    'scipy',
    'matplotlib==3.4.3',     
    'pandas==1.20.3',
    'f90wrap'
]

for pkg in python_reqs:
    os.system('pip install ' + pkg)

# Cloning repositories
os.system('mkdir ../dmk_utilities')
repos_list = ['dmk_solver', 'globals', 'linear_algebra', 'p1galerkin', 'geometry']
for repo in repos_list:
    os.system('cd ../dmk_utilities && git clone https://gitlab.com/enrico_facca/' + repo + '.git')

## Creating builds
os.system('cd ../dmk_utilities/dmk_solver && mkdir build/ && cd build && cmake ../ && make')
os.system('cd ../dmk_utilities/dmk_solver/build && make f2py_interface_dmk')

# Getting repository's location
cwd_ = os.getcwd()
f = open("nextrout_location.txt", "w")
f.write(cwd_)
f.close()

print('Installation finished!')
print('In case of errors, please check the Troubleshooting section at https://gitlab.com/enrico_facca/dmk_solver')
