#!/bin/bash

git config --global credential.helper store

for repo2 in python_scripts
do
	link='https://github.com/Danielaleite/'${repo2}'.git'
	if [ ! -d "${repo2}" ]; then
		echo "Cloning $THEPATH ( $repo2)"
		git clone $link
	else
		echo "Pulling $repo2"
		cd "$repo2" && git checkout master && git pull && cd ~-
	fi

done