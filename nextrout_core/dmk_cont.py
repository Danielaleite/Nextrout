### import needed stuff
# 
# Loading stardard and dmk pyhton modules
#

# Standard tools
import sys
import numpy as np
from itertools import combinations


# Import I/O for timedata
try:
    sys.path.append('../dmk_utilities/globals/python/timedata/')
    import timedata as td
except:
    print("Global repo non found")

# Import geometry tools
sys.path.append('../dmk_utilities/geometry/python/')
import meshtools as mt
sys.path.append('../dmk_utilities/dmk_solver/otp_solver/preprocess/assembly/')
import example_grid

# Import dmk tools
sys.path.append('../dmk_utilities/dmk_solver/otp_solver/python/')
import dmk_p1p0 
sys.path.append('/../dmk_utilities/dmk_solver/build/python/fortran_python_interface/')
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

def grid_gen(ndiv, nref=1, flag_grid='unitsquare'):
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
    centers= {}
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

    if 'dirac' in forcing_flag:
        
        Nplus = extra_info['Nplus']
        Nminus = extra_info['Nminus']

        fplus = extra_info['fplus']
        fminus=  extra_info['fminus']

        xplus = extra_info['xplus'] # [[0.1,0.2],[0.3,0.4],[0.1,0.7]]
        xminus = extra_info['xminus'] # [[0.6,0.2],[0.8,0.4]]

        # set array forcing_dirac "evoluation" f=f^{+}-f^{-} on grid nodes
        forcing=np.zeros(grid.nnode)
        for i in range(Nplus):
            inode=mt.Inode(coord,xplus[i])
            forcing[inode]=fplus[i]
        for i in range(Nminus):
            inode=mt.Inode(coord,xminus[i])
            forcing[inode]=-fminus[i]

    elif 'rect' in forcing_flag:

        rectangles_source = extra_info[0] # it should be a list of lists: every sublist has 3 elements: (x,y), w, h
        rectangles_sink = extra_info[1]
        
        forcing =np.zeros(grid.nnode)
        # iterate over the center of the triangles
        for cent in centers.keys():

          x,y = centers[cent]

          for rect in rectangles_source:

            x_source = rect[0][0]
            y_source = rect[0][1]
            wo = rect[1]
            ho = rect[2]

            if (x_source <= x <= x_source+wo) and (y_source <= y <= y_source+ho):
              inode=mt.Inode(coord,centers[cent])
              #print(inode)
              forcing[inode]=1

          for rect in rectangles_sink:

            x_sink = rect[0][0]
            y_sink = rect[0][1]
            wi = rect[1]
            hi = rect[2]

            if (x_sink <= x <= x_sink+wi) and (y_sink <= y <= y_sink+hi):
              inode=mt.Inode(coord,centers[cent])
              #print(inode)
              forcing[inode]=-1
    else:
        raise ValueError('forcing not defined')

    return forcing

def kappa_generator(coord, kappa_flag):

    if kappa_flag == 'random':
        x=coord[0]; y=coord[1]
        fvalue=0.1+np.random.random_sample()
    else:
        raise ValueError('kappa not defined!')

    return fvalue


def dmk_cont(forcing, beta_c, ndiv, tdens0 = None, nref= 0, flag_grid = 'unitsquare', kappa_flag = None):

    length=1.0/float(ndiv)

    # build grid using prebuild examples 
    points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)

    # initialized fortran variable for the spatial discretization
    [grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
    ncell=grid.ncell
    ntdens=grid.ncell
    npot=subgrid.nnode

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
    Dmkcontrols.get_from_file(ctrl,'dmk.ctrl')
    ctrl.fn_tdens='tdens.dat'
    ctrl.fn_pot='pot.dat'
    ctrl.fn_statistics='dmk.log'
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

    return tdpot

