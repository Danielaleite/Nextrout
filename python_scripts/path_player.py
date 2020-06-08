import sys
import os

# let's try to access to the different folders and files without adding relative paths and special dependencies, just
# abs paths.

def get_root_path(path):
    current_path = os.path.abspath(path)
    path_parts = current_path.split('/')
    root = '/'
    for word in path_parts:
        if word != 'Nextrout':
            root =root+word+'/'
        else:
            root=root+word+'/'
            break
    return root

#get_root_path('Nextrout')

print(os.path.abspath('/Nextrout'))