### import needed stuff
# 
# Loading stardard and dmk pyhton modules
#

# Standard tools
import sys
from pathlib import Path

import numpy as np
from itertools import combinations
import os

# Accessing root path
file_path = Path(__file__).resolve().parent
root = file_path.parents[1]
print(root)
# Import I/O for timedata
try:
    sys.path.append(str(root / 'dmk_utilities/globals/python/timedata/'))
    import timedata as td
except:
    print("Global repo non found")

# Import geometry tools
sys.path.append(str(root / 'dmk_utilities/geometry/python/'))
import meshtools as mt
sys.path.append(str(root / 'dmk_utilities/dmk_solver/otp_solver/preprocess/assembly/'))
import example_grid

# Import dmk tools
sys.path.append(str(root / 'dmk_utilities/dmk_solver/otp_solver/python/'))
import dmk_p1p0 
sys.path.append(str(root / 'dmk_utilities/dmk_solver/build/python/fortran_python_interface/'))
from dmk import (Dmkcontrols,    # controls for dmk simulations)
                 Timefunctionals, # information of time/algorithm evolution
                Dmkinputsdata, # structure variable containg inputs data
                 build_subgrid_rhs, #procedure to preprocess forcing term f
                 Tdenspotentialsystem, # input/output result tdens, pot
                dmkp1p0_steady_data   # main interface subroutine
                )
# Import plot tools
import matplotlib.pyplot as plt
import matplotlib.tri as mtri

def grid_gen(ndiv, nref=1, flag_grid='rect_cnst'):
    #
    # Define mesh for spatial disctetization.
    # Build the "coord" and "topol" numpy arrays describing coordinate and topology of the mesh.
    #

    # set mesh size 
    length=1.0/float(ndiv)

    # build grid using prebuild examples 
    points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)

    # initialized fortran variable for the spatial discretization
    [grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)

    return grid, subgrid, points, vertices, coord,topol,element_attributes

def center_computation(topol, coord):
    centers= np.zeros((len(topol),2))
    coordinates = coordinate_computation(coord)
    k=-1
    for T in topol: #T is a triangle
      k+=1
      edges_in_T = list(combinations(T, 2))
      coord_of_nodes_in_T = np.array([coordinates[str(node)] for node in T])
      center = sum(coord_of_nodes_in_T)/3
      centers[k] = center
    return centers

def coordinate_computation(coord):
    coordinates = {}
    k=-1
    for index_ in coord:
      k+=1
      coordinates [str(k)] = index_[:2]
    return coordinates

def forcing_generator(forcing_flag, grid, coord, topol,extra_info):

    centers = center_computation(topol, coord)

    grid_source_indices = []
    grid_sink_indices = []

    if 'dirac' in forcing_flag:
        
        Nplus = extra_info['Nplus']
        Nminus = extra_info['Nminus']

        fplus = extra_info['fplus']
        fminus=  extra_info['fminus']

        xplus = extra_info['xplus'] # [[0.1,0.2],[0.3,0.4],[0.1,0.7]]
        xminus = extra_info['xminus'] # [[0.6,0.2],[0.8,0.4]]

        # set array forcing_dirac "evolution" f=f^{+}-f^{-} on grid nodes
        forcing=np.zeros(grid.nnode)
        for i in range(Nplus):
            inode=mt.Inode(coord,xplus[i])
            forcing[inode]=1#fplus[i]
            grid_source_indices.append(inode)
        for i in range(Nminus):
            inode=mt.Inode(coord,xminus[i])
            forcing[inode]=-1#fminus[i]
            grid_sink_indices.append(inode)

    elif 'rect' in forcing_flag:

        rectangles_source = extra_info[0] # it should be a list of lists: every sublist has 3 elements: (x,y), w, h
        rectangles_sink = extra_info[1]
        
        forcing =np.zeros(grid.nnode)
        # iterate over the center of the triangles
        for cent in centers:

          x,y = cent

          for rect in rectangles_source:

            x_source = rect[0][0]
            y_source = rect[0][1]
            wo = rect[1]
            ho = rect[2]

            if (x_source <= x <= x_source+wo) and (y_source <= y <= y_source+ho):
              
              inode = mt.Inode(centers,cent)
              #print(inode)
              forcing[inode]=1
              grid_source_indices.append(inode)

          for rect in rectangles_sink:

            x_sink = rect[0][0]
            y_sink = rect[0][1]
            wi = rect[1]
            hi = rect[2]

            if (x_sink <= x <= x_sink+wi) and (y_sink <= y <= y_sink+hi):
              
              inode = mt.Inode(centers,cent)
              #print(inode)
              forcing[inode]=-1
              grid_sink_indices.append(inode)

    else:
        raise ValueError('forcing not defined')

    ### triangle indices
    triang_source_indices = []
    triang_sink_indices = []
    k=-1
    for T in topol:
        k+=1
        if T[0] in grid_source_indices or T[1] in grid_source_indices or T[2] in grid_source_indices:
            triang_source_indices.append(k)
        elif T[0] in grid_sink_indices or T[1] in grid_sink_indices or T[2] in grid_sink_indices:
            triang_sink_indices.append(k)

    npos = len(forcing[forcing>0])
    nneg = len(forcing[forcing<0])

    for i in range(len(forcing)):
        if forcing[i]>0:
            forcing[i] = forcing[i]/npos
        elif forcing[i]<0:
            forcing[i] = forcing[i]/nneg

    assert sum(forcing)<0.01
    
    return forcing, triang_source_indices,triang_sink_indices

def kappa_generator(coord, kappa_flag):

    if kappa_flag == 'random':
        x=coord[0]; y=coord[1]
        fvalue=0.1+np.random.random_sample()
    else:
        raise ValueError('kappa not defined!')

    return fvalue


def dmk_cont(forcing, beta_c, ndiv, extra_info, niter=80,tdens0 = None, nref= 0, flag_grid =
'rect_cnst', kappa_flag = None, storing = None):

    if storing == None:
        storing = '.'

    length=1.0/float(ndiv)

    # build grid using prebuild examples 
    points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)

    # initialized fortran variable for the spatial discretization
    [grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
    ncell=grid.ncell
    ntdens=grid.ncell
    npot=subgrid.nnode
    
    # generate forcing

    # forcing, triang_source_indices, triang_sink_indices = forcing_generator(
    #     forcing, grid, coord, topol, extra_info=extra_info
    # )

    # initial integrated forcing term
    rhs=np.zeros(subgrid.ncell)
    dmk_p1p0.build_subgrid_rhs(subgrid,rhs, np.zeros(ncell),forcing)

    # Init and set "container" with inputs for dmk simulation
    dmkin=Dmkinputsdata.DmkInputs()
    Dmkinputsdata.dmkinputs_constructor(dmkin,0,ntdens,npot,True) # this True set to default all varaibles

    # integrate forcing term w.r.t. p1 base function
    build_subgrid_rhs(subgrid, dmkin.rhs, np.zeros(grid.ncell),forcing)
    dmkin.pflux = beta_c

    # Init "container" variable with tdens(mu) and potential(u) varaible
    tdpot=Tdenspotentialsystem.tdpotsys()
    Tdenspotentialsystem.tdpotsys_constructor(tdpot,0,ntdens, npot,1)
    tdpot.tdens[:]=1.0
    ### defining tdens0

    if tdens0 is None:
        tdpot.tdens[:]=1.0
    else:
        tdpot.tdens = tdens0

    ### defining kappa
    if kappa_flag is not None:
        # compute functions on cell centroids
        ncell=len(topol)
        bar_cell=mt.make_bar(coord,topol).transpose()
        kappa_cell=np.zeros([ncell]);
        for i in range(ncell):
            kappa_cell[i] = kappa_generator(bar_cell[:,i])

    # init and set controls
    ctrl = Dmkcontrols.DmkCtrl()
    Dmkcontrols.get_from_file(ctrl,root / 'Nextrout/nextrout_core/dmk_cont.ctrl')
    ctrl.fn_tdens=storing+'/tdens.dat'
    ctrl.fn_pot=storing+'/pot.dat'
    ctrl.fn_statistics=storing+'/dmk.log'
    ctrl.max_time_iterations = niter
    #
    # init type for storing evolution/algorithm info
    #
    timefun=Timefunctionals.evolfun()
    Timefunctionals.evolfun_constructor(timefun, 0,
                                            ctrl.max_time_iterations,
                                            ctrl.max_nonlinear_iterations)

    # solve with dmk
    info=0
    id_subgrid=np.ones(1,np.int32)
    dmkp1p0_steady_data(grid, subgrid,  id_subgrid[0],tdpot, dmkin, ctrl, info,timefun=timefun)

    if info ==0: print('convergence achieved!.')

    return tdpot, timefun



'''
ndiv = 20

length=1.0/float(ndiv)

# set grid example
flag_grid='unitsquare'

# build grid using prebuild examples 
points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)

# initialized fortran variable for the spatial discretization
[grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
ncell=grid.ncell
ntdens=grid.ncell
npot=subgrid.nnode


Nplus=3
Nminus=2

fplus=[1,2,3]
fminus=[4,2]

xplus=[[0.1,0.2],[0.3,0.4],[0.1,0.7]]
xminus=[[0.6,0.2],[0.8,0.4]]

# set array forcing_dirac "evoluation" f=f^{+}-f^{-} on grid nodes
forcing=np.zeros(grid.nnode)
for i in range(Nplus):
    inode=mt.Inode(coord,xplus[i])
    forcing[inode]=fplus[i]
for i in range(Nminus):
    inode=mt.Inode(coord,xminus[i])
    forcing[inode]=-fminus[i]

# initial integrated forcing term
rhs=np.zeros(subgrid.ncell)
dmk_p1p0.build_subgrid_rhs(subgrid,rhs, np.zeros(ncell),forcing)

# Init and set "container" with inputs for dmk simulation
dmkin=Dmkinputsdata.DmkInputs()
Dmkinputsdata.dmkinputs_constructor(dmkin,0,ntdens,npot,True) # this True set to default all varaibles

# integrate forcing term w.r.t. p1 base function
build_subgrid_rhs(subgrid, dmkin.rhs, np.zeros(grid.ncell),forcing)
dmkin.pflux = 1#beta_c

# Init "container" variable with tdens(mu) and potential(u) varaible
tdpot=Tdenspotentialsystem.tdpotsys()
Tdenspotentialsystem.tdpotsys_constructor(tdpot,0,ntdens, npot,1)
tdpot.tdens[:]=1.0
### defining tdens0
tdens0 = None
if tdens0 is None:
    tdpot.tdens[:]=1.0
else:
    tdpot.tdens = tdens0

### defining kappa
kappa_flag = None
if kappa_flag is not None:
    # compute functions on cell centroids
    ncell=len(topol)
    bar_cell=mt.make_bar(coord,topol).transpose()
    kappa_cell=np.zeros([ncell]);
    for i in range(ncell):
        kappa_cell[i] = kappa_generator(bar_cell[:,i])

# init and set controls
ctrl = Dmkcontrols.DmkCtrl()
Dmkcontrols.get_from_file(ctrl,root+'/nextrout_core/dmk_cont.ctrl')
ctrl.fn_tdens=storing+'/tdens.dat'
ctrl.fn_pot=storing+'/pot.dat'
ctrl.fn_statistics=storing+'/dmk.log'
ctrl.max_time_iterations = 300
#
# init type for storing evolution/algorithm info
#
timefun=Timefunctionals.evolfun()
Timefunctionals.evolfun_constructor(timefun, 0,
                                        ctrl.max_time_iterations,
                                        ctrl.max_nonlinear_iterations)

# solve with dmk
info=0
dmkp1p0_steady_data(grid, subgrid, tdpot, dmkin, ctrl, info,timefun=timefun)

if info ==0: print('convergence achieved!.')


'''
