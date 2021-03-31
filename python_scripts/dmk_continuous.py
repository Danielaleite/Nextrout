# 
# Loading stardard and dmk pyhton modules
#

# Standard tools
import sys
import numpy as np
!pip install pymesh
!pip install numpy


# Import I/O for timedata
try:
    sys.path.append('globals/python/timedata/')
    import timedata as td
except:
    print("Global repo non found")

# Import geometry tools
sys.path.append('geometry/python/')
import meshtools as mt
sys.path.append('dmk_solver/otp_solver/preprocess/assembly/')
import example_grid

# Import dmk tools
sys.path.append('dmk_solver/otp_solver/python')
import dmk_p1p0 
sys.path.append('dmk_solver/build/python/fortran_python_interface/')
from dmk import (Dmkcontrols,    # controls for dmk simulations)
                 Timefunctionals, # information of time/algorithm evolution
                Dmkinputsdata, # structure variable containg inputs data
                 build_subgrid_rhs, #procedure to preprocess forcing term f
                 Tdenspotentialsystem, # input/output result tdens, pot
                dmkp1p0_steady_data   # main interface subroutine
                )
# Import plot tools
import matplotlib.pyplot as plt

sys.path.append('geometry/python/')
import meshtools as mt

def initializing(ndiv,nref, flag_grid = 'unitsquare'):

	# set mesh size 
	ndiv=18
	length=1.0/float(ndiv)
	nref=0

	# set grid example
	flag_grid='unitsquare'

	# build grid using prebuild examples 
	points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)

	# initialized fortran variable for the spatial discretization
	[grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
	ncell=grid.ncell
	ntdens=grid.ncell
	npot=subgrid.nnode

	#
	# set number, value and location of source and sink points
	# 
	Nplus=3
	Nminus=2

	fplus=[1,2,3]
	fminus=[4,2]

	xplus=[[0.1,0.2],[0.3,0.4],[0.1,0.7]]
	xminus=[[0.6,0.2],[0.8,0.4]]

	# set array forcing_dirac "evoluation" f=f^{+}-f^{-} on grid nodes
	forcing_dirac=np.zeros(grid.nnode)
	for i in range(Nplus):
	    inode=mt.Inode(coord,xplus[i])
	    forcing_dirac[inode]=fplus[i]
	for i in range(Nminus):
	    inode=mt.Inode(coord,xminus[i])
	    forcing_dirac[inode]=-fminus[i]

	# plot location oand intensity of source and sink
	ax = plt.gca()
	ax.cla() 
	ax.set_xlim((0, 1))
	ax.set_ylim((0, 1))
	# some data
	for i in range(Nplus):
	    print((xplus[i][0],xplus[i][1]))
	    circle = plt.Circle((xplus[i][0],xplus[i][1]),0.01*fplus[i], color='red', fill=True)
	    ax.add_patch(circle)
	for i in range(Nminus):
	    circle = plt.Circle((xminus[i][0],xminus[i][1]),0.01*fminus[i], color='blue', fill=True)
	    ax.add_patch(circle)
	plt.show()

	    

	# initial integrated forcing term
	rhs=np.zeros(subgrid.ncell)
	dmk_p1p0.build_subgrid_rhs(subgrid,rhs, np.zeros(ncell),forcing_dirac)


### before doing this, config the dmk utilities!
print()