#!/bin/bash

git config --global credential.helper store

for repo in globals linear_algebra geometry vtk p1galerkin muffe_p1p0 muffe_sparse_optimization spectral
do 
    link='https://gitlab.com/enrico_facca/'${repo}'.git'
    if [ ! -d "${repo}" ]; then
        echo "Cloning $THEPATH ( $repo )"
        git clone $link
    else
        echo "Pulling $repo"
        cd "$repo" && git checkout master && git pull && cd ~-
    fi
done
    
    