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

def grid_gen(ndiv, nref=0, flag_grid='unitsquare'):
    #
    # Define mesh for spatial disctetization.
    # Build the "coord" and "topol" numpy arrays describing coordinate and topology of the mesh.
    #

    # set mesh size 
    length=1.0/float(ndiv)
    nref=0

    # build grid using prebuild examples 
    points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)

    # initialized fortran variable for the spatial discretization
    [grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
    ncell=grid.ncell
    ntdens=grid.ncell
    npot=subgrid.nnode

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

    if forcing_flag == 'dirac':
        
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

    elif forcing_flag == 'rect_cnst':

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


def dmk_cont(ndiv, forcing, beta, tdens0 = None,  nref= 0, flag_grid = 'unitsquare', kappa_flag = None):

    ### building the grid
    # set mesh size 
    length=1.0/float(ndiv)

    # build grid using prebuild examples 
    points, vertices, coord,topol,element_attributes = example_grid.example_grid(flag_grid,length)
    # initialized fortran variable for the spatial discretization
    [grid,subgrid]=dmk_p1p0.init_geometry(topol, coord, 1)
    ncell=grid.ncell
    ntdens=grid.ncell
    npot=subgrid.nnode

    ### defining dmk controls
    ctrl = Dmkcontrols.DmkCtrl()
    ## redef param.
    ctrl.tolerance_system_variation = .000001
    print(ctrl.tolerance_system_variation)

    # Init and set "container" with inputs for dmk simulation
    dmkin=Dmkinputsdata.DmkInputs()
    Dmkinputsdata.dmkinputs_constructor(dmkin,0,ntdens,npot,True) # this True set to default all varaibles

    ## redef param.
    dmkin.pflux = beta
    print(dmkin.pflux)

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

    ### defining rhs

    # initial integrated forcing term
    rhs=np.zeros(subgrid.ncell)
    # integrate forcing term w.r.t. p1 base function
    build_subgrid_rhs(subgrid, dmkin.rhs, np.zeros(grid.ncell),forcing)

    # init and set controls
    ctrl = Dmkcontrols.DmkCtrl()
    Dmkcontrols.get_from_file(ctrl,'dmk.ctrl')
    ctrl.fn_tdens='tdens.dat'
    ctrl.fn_pot='pot.dat'
    ctrl.fn_statistics='dmk.log'

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

#generate grid
ndiv = 25
grid, subgrid, points, vertices, coord,topol,element_attributes = grid_gen(ndiv)

# generate forcing

forcing_flag = 'rect_cnst'

if forcing_flag == 'rect_cnst':

    x_source1, y_source1 = (.2,.2)
    x_source2, y_source2 = (.2,.7)
    wo1 = .05
    ho1 = .1
    rectangles_source = [[(x_source1,y_source1),wo1,ho1],[(x_source2,y_source2),wo1,ho1]] # bottom left cornner, width, height

    x_sink, y_sink = (.8,.8)
    wi = .1
    hi = .1
    rectangles_sink = [[(x_sink,y_sink),wi,hi]]

    extra_info = [rectangles_source,rectangles_sink]

    

elif forcing_flag == 'dirac':
    Nplus=3
    Nminus=2

    fplus=[1,2,3]
    fminus=[4,2]

    xplus=[[0.1,0.21],[0.3,0.4],[0.1,0.7]]
    xminus=[[0.6,0.2],[0.8,0.4]]

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}


forcing = forcing_generator(forcing_flag, grid, coord, topol, extra_info=extra_info)
triang = mtri.Triangulation(coord.transpose()[0,:], coord.transpose()[1,:], topol)

fig, ax = plt.subplots()
len(coord.transpose()[0,:])
ax.tricontour(triang, forcing, levels=40, linewidths=0.1, cmap = 'RdBu_r')
#cntr2 = ax.tricontourf(coord.transpose()[0,:], coord.transpose()[1,:], forcing_dirac, levels=14, cmap="RdBu_r")

plt.subplots_adjust(hspace=0.5)
plt.show()

#run dmk
beta = 1.5
tdpot = dmk_cont(ndiv,forcing,beta)

# plotting results
#triang = mtri.Triangulation(coord.transpose()[0,:], coord.transpose()[1,:], topol)
fig1, ax1 = plt.subplots(figsize=(8, 8)); ax1.set_aspect('equal')
tpc = ax1.tripcolor(triang, tdpot.tdens , cmap='Greens')
fig1.colorbar(tpc)
ax1.set_title('Optimal transportation density $\mu^*$')
plt.show()