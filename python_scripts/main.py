# -*- coding: utf-8 -*-
#!/usr/bin/env python

#######################################################
# This program takes 
# f_ctrl        : inputs file wiht flag/controls/ path 
#                 See 'inputs.ctrl' in ../script_controls/
# dir_inputs    : folder where store all inputs data
# dir_inputs    : folder where store data in vtk format
#                
#######################################################

import numpy as np
import subprocess
import sys
import os

#
# geometry
#
#sys.path.insert(0,'./geometry2d')
import meshtools as mt
import example_grid as ex_grid
from pyvtk import *



#
# inputs
#

import example_forcing as ex_forcing
import make_dirichlet as dirichlet
import make_optimal_tdens as make_optdens

import common as common
import timecell as timecell
import cell as cell
import scalars as scalars




nnodeincell = 3
ndim = 2


############################################
# read control path
f_ctrl=sys.argv[1]
dir_inputs=sys.argv[2]
dir_inputs_vtk=sys.argv[3]


if ( dir_inputs == ''):
    dir_inputs='./'
if ( dir_inputs_vtk == ''):
    dir_inputs_vtk='./'

# set grid destination
file_grid=dir_inputs+'/'+'grid.dat'
file_grid_vtk=dir_inputs_vtk+'/'+'grid'
file_poly_vtk=dir_inputs_vtk+'/'+'poly'

# set inputs destinations
filename_forcing=dir_inputs+'/'+'forcing.dat'
filename_source=dir_inputs+'/'+'source.dat'
filename_sink=dir_inputs+'/'+'sink.dat'
filename_dirac_forcing=dir_inputs+'/'+'dirac_forcing.dat'
filename_dirac_source=dir_inputs+'/'+'dirac_source.dat'
filename_dirac_sink=dir_inputs+'/'+'dirac_sink.dat'


# cell time data (changing in time)
filename_kappa=dir_inputs+'/'+'kappa.dat'

# scalar time data (changing in time)
filename_pflux=dir_inputs+'/'+'pflux.dat'
filename_pmass=dir_inputs+'/'+'pmass.dat'
filename_decay=dir_inputs+'/'+'decay.dat'


# cell data (fix in time)
filename_measure_source=dir_inputs+'/'+'measure_source.dat'
filename_measure_sink=dir_inputs+'/'+'measure_sink.dat'
filename_tdens0=dir_inputs+'/'+'tdens0.dat'
filename_optdens=dir_inputs+'/optdens.dat'

#vtk
filename_inputs_vtk=dir_inputs_vtk+'/'+'inputs'


############################################
# reads flags and other values
ctrl = open(f_ctrl, "r")
ctrl_lines = ctrl.readlines()

# flags grid
flag_grid=str(common.read_column('flag_grid',0,ctrl_lines))
extra_grid=str(common.read_column('flag_grid',1,ctrl_lines))
ndiv=int(common.read_column('ndiv',0,ctrl_lines))
nref=int(common.read_column('nref',0,ctrl_lines))

# flag_inputs
flag_source=str(common.read_column('flag_source',0,ctrl_lines))
extra_source=str(common.read_column('flag_source',1,ctrl_lines))

flag_sink=str(common.read_column('flag_sink',0,ctrl_lines))
extra_sink=str(common.read_column('flag_sink',1,ctrl_lines))

flag_normalize=int(common.read_column('flag_normalize',0,ctrl_lines))

flag_dirichlet=str(common.read_column('flag_dirichlet',0,ctrl_lines))
extra_dirichlet=str(common.read_column('flag_dirichlet',1,ctrl_lines))

flag_neumann=str(common.read_column('flag_neumann',0,ctrl_lines))
extra_neumann=str(common.read_column('flag_neumann',1,ctrl_lines))

# grid dependent quantities
flag_tdens0=str(common.read_column('flag_tdens0',0,ctrl_lines))
extra_tdens0=str(common.read_column('flag_tdens0',1,ctrl_lines))

flag_kappa=str(common.read_column('flag_kappa',0,ctrl_lines))
extra_kappa=str(common.read_column('flag_kappa',1,ctrl_lines))

# scalar quantities
flag_pflux=str(common.read_column('flag_pflux',0,ctrl_lines))
extra_pflux=str(common.read_column('flag_pflux',1,ctrl_lines))

flag_pmass=str(common.read_column('flag_pmass',0,ctrl_lines))
extra_pmass=str(common.read_column('flag_pmass',1,ctrl_lines))

flag_decay=str(common.read_column('flag_decay',0,ctrl_lines))
extra_decay=str(common.read_column('flag_decay',1,ctrl_lines))


# time interval
tzero=float(common.read_column('tzero',0,ctrl_lines))
tmax=float(common.read_column('tmax',0,ctrl_lines))



sources=[];
sinks=[];
dirac_sources=[];
dirac_sinks=[];
kappas=[];
tdens0_functions=[];
optdens_functions=[];
dirac_points=[];


length=0.0

###############################
# build mesh
if ( flag_grid != 'use'):
    length=1.0/float(ndiv)
    points, vertices, coordinates,elements,element_attributes = ex_grid.example_grid(flag_grid,length,extra_grid)    
    
    coord=np.array(coordinates)
    print (coord.shape)
    if ( coord.shape[1] == 2 ):
        zcoord = np.zeros([coord.shape[0],1])
        coord=np.append(coord, zcoord, axis=1)
    
    topol=np.array(elements)
    #print (element_attributes)
    try:
        print (len(element_attributes))
        for i in range(len(element_attributes)):
            print (i, elements_attributes[i])

        flags=np.array(element_attributes)
    except NameError:
        flags=np.int_(range(len(topol)))
    

    
    mt.write_poly(points,vertices,str(file_poly_vtk),'vtk')
else:
    # reads new coord topol and flags of refined grid
    coord, topol, flags = mt.read_grid(extra_grid)
    if ( coord.shape[1] == 2 ):
        zcoord = np.zeros([coord.shape[0],1])
        coord=np.append(coord, zcoord, axis=1)

    


#
# call external program to refine mesh
#
if (nref > 0):
    file_grid0=dir_inputs+'/'+'grid0.dat'
    file_grid0_vtk=dir_inputs_vtk+'/'+'grid0'
    mt.write_grid(coord,topol,file_grid0,'dat',flags)
    mt.write_grid(coord,topol,file_grid0_vtk,'vtk',flags)
    cmd='./uniform_refinement/code/uniform_refinement.out'+\
        ' '+str(file_grid0)+' '+str(file_grid)+' '+str(nref)
    subprocess.call(cmd,shell=True)
    # reads new coord topol and flags of refined grid
    coord, topol, flags = mt.read_grid(file_grid)
    if ( coord.shape[1] == 2 ):
        zcoord = np.zeros((coord.shape[0],1))
        coord = np.append(coord, zcoord, axis=1)


mt.write_grid(coord,topol,file_grid,'dat',flags)
mt.write_grid(coord,topol,str(file_grid_vtk),'vtk',flags)
    
# computing basic geomtry info of the mesh
ncell=len(topol)
nnode=len(coord)
size_cell=mt.make_size(coord,topol)
bar_cell=mt.make_bar(coord,topol)

if ( length == 0.0 ):
    length = np.sqrt( 2.0* np.sum(size_cell) / float( len(size_cell) ) )
else:
    length = length / float(nref+1)

ncell = len(topol)
nnode = len(coord)


###############################
# prepare sources ans sink functions if not defined
# append sources or dirac nodes
if ( ( len(sources) == 0) and (len(dirac_sources) == 0) ) :
    ex_forcing.example_source(sources,dirac_sources,
                              str(flag_source),str(extra_source))

# append sinks or driac nodes
if ( ( len(sinks) == 0) and (len(dirac_sinks) == 0) ) :
    ex_forcing.example_sink(sinks,dirac_sinks,
                        str(flag_sink),str(extra_sink))


##########################
# forcing init.
source_tria=np.zeros([ncell,1])
sink_tria=np.zeros([ncell,1])
forcing_tria=np.zeros([ncell,1])
dirac_forcing_nodes=np.zeros([nnode,1])
dirac_source_nodes=np.zeros([nnode,1])
dirac_sink_nodes=np.zeros([nnode,1])


steady_source=False
steady_sink=False
steady_forcing=False
dirac_source_steady=False
dirac_sink_steady=False
dirac_forcing_steady=False

if ( True ) :
    time=tzero
   

    # open file and write dimensions
    file_forcing=open(filename_forcing, 'w')
    file_source=open(filename_source, 'w')
    file_sink=open(filename_sink, 'w')

    timecell.write2file(file_forcing,time,True,
                        forcing_tria)
    timecell.write2file(file_source,time,True,
                        forcing_tria)
    timecell.write2file(file_sink,time,True,
                        forcing_tria)

    # open file and write dimensions
    file_dirac_forcing=open(filename_dirac_forcing, 'w')
    file_dirac_source=open(filename_dirac_source, 'w')
    file_dirac_sink=open(filename_dirac_sink, 'w')

    timecell.write2file(file_dirac_forcing,time,True,
                        dirac_forcing_nodes)
    timecell.write2file(file_dirac_source,time,True,
                        dirac_source_nodes)
    timecell.write2file(file_dirac_sink,time,True,
                        dirac_sink_nodes)
        
    #cycle in time
    while ( (time <= tmax) and 
            ( (not steady_forcing ) or (not steady_dirac_forcing) )
        ):
        #make source
        source_tria, steady_source=ex_forcing.make_source(
            source_tria,steady_source,
            sources,
            time,flags,bar_cell)
    

        dirac_source_nodes, dirac_source_steady= ex_forcing.make_dirac_source(
            dirac_source_nodes,dirac_source_steady,dirac_sources,time,
            str(extra_source),coord)
        
        if ( str(flag_source) == 'use') :
            cmd='cp '+ str(extra_source) +' '+str(filename_source)
            subprocess.call(cmd,shell=True)
            steady_source=True


        #make sink
        sink_tria, steady_sink=ex_forcing.make_sink(
            sink_tria,steady_sink,
            sinks,
            time,flags,bar_cell)

        dirac_sink_nodes, dirac_sink_steady = ex_forcing.make_dirac_sink(
            dirac_sink_nodes,dirac_sink_steady,dirac_sinks,time,
            str(extra_source),coord)

        if ( str(flag_source) == 'use') :
            cmd='cp '+ str(extra_source) +' '+str(filename_source)
            subprocess.call(cmd,shell=True)
            steady_sink=True


        steady_forcing=(
            steady_source and 
            steady_sink)
        steady_dirac_forcing=( 
            dirac_source_steady and 
            dirac_sink_steady)


        if ( flag_dirichlet =='no') :
            print( 'correction' )
            source_tria,sink_tria, dirac_source_nodes,dirac_sink_nodes=ex_forcing.correction(
                topol,
                source_tria,sink_tria,
                dirac_source_nodes,dirac_sink_nodes,
                size_cell)

        if (flag_normalize == 1):
            source_tria,sink_tria,dirac_source_nodes,dirac_sink_nodes=ex_forcing.normalize( source_tria,sink_tria, dirac_source_nodes, dirac_sink_nodes,size_cell)
            

        forcing_tria = source_tria-sink_tria
        dirac_forcing_nodes = dirac_source_nodes-dirac_sink_nodes
    
        timecell.write2file( 
            file_forcing,time,False,
            forcing_tria,steady_forcing)
        timecell.write2file(
            file_source,time,False,
            source_tria,steady_source)
        timecell.write2file(
            file_sink,time,False,
            sink_tria,steady_sink)

        timecell.write2file(
            file_dirac_forcing,time,False,
            dirac_forcing_nodes,steady_dirac_forcing)
        timecell.write2file(
            file_dirac_source,time,False,
            dirac_source_nodes,dirac_source_steady)
        timecell.write2file(
            file_dirac_sink,time,False,
            dirac_sink_nodes,dirac_sink_steady)
    
        time=common.next_time(time)
        
        ### close files
        file_forcing.close()
        file_source.close()
        file_sink.close()
        file_dirac_forcing.close()
        file_dirac_source.close()
        file_dirac_sink.close()

        #####################
        # write last measure
        out=np.zeros([ncell,4])
        out[:,0:3]=bar_cell[:,:]
        
        out[:,3]=source_tria[:][0]*size_cell[:]
        cell.write2file(filename_measure_source,out)
    
        out[:,3]=sink_tria[:][0]*size_cell[:]
        cell.write2file(filename_measure_sink,out)
        
if ( str(flag_source) == 'use') :
    cmd='cp '+ str(extra_source) +' '+str(filename_source)
    subprocess.call(cmd,shell=True)

if ( str(flag_sink) == 'use') :    
    cmd='cp '+ str(extra_sink) +' '+str(filename_sink)
    subprocess.call(cmd,shell=True)







# #################################################################
# kappa
if ( not kappas):
    kappas=[]
    kappas.append(timecell.example(flag_kappa,extra_kappa))

# definitions 
kappa_tria=np.zeros([ncell,1])
steady=False
time=tzero

# open file and write dimensions
file_out=open(filename_kappa, 'w')
timecell.write2file(file_out,time,True,kappa_tria)

#cycle in time
while (time <= tmax) and (not steady ):
    # evalute 
    kappa_tria, steady = timecell.build(
        kappa_tria, steady, 
        kappas,time,bar_cell,flags)
    # write 2 file
    timecell.write2file(file_out,time,False,kappa_tria,steady)
    
    # next time
    time=common.next_time(time)
file_out.close()

#####################################################################
# tdens0
#####################################################################
if ( not tdens0_functions):
    tdens0_functions=[]
    tdens0_functions.append(timecell.example(flag_tdens0,extra_tdens0))

# definitions 
tdens0_tria=np.zeros([ncell,1])
steady=False
time=tzero

# open file and write dimensions
file_out=open(filename_tdens0, 'w')
timecell.write2file(file_out,time,True,tdens0_tria)

#cycle in time
while (time <= tmax) and (not steady ):
    # evalute 
    tdens0_tria, steady = timecell.build(
        tdens0_tria, steady, 
        tdens0_functions,time,bar_cell,flags)
    # write 2 file
    timecell.write2file(file_out,time,False,tdens0_tria,steady)
    
    # next time
    time=common.next_time(time)
file_out.close()


#####################################################################
# Dirichlet 
#####################################################################
# if ( flag_dirichlet != 'no'):
#     filename=dir_inputs+'/dirichlet.dat'
#     dirichlet.example_dirichlet(flag_dirichlet,coord,extra_dirichlet,filename)
    


#####################################################################
# OPTDENS optimal transpost density
#####################################################################
if ( not optdens_functions):
    optdens_functions=[]
if (len(optdens_functions) == 0):
    optdens_functions.append(make_optdens.define_optimal_tdens(
        flag_grid, flag_source, flag_sink, 
        extra_grid, extra_source, extra_sink ))

try:
    # definitions 
    optdens_tria=np.zeros([ncell,1])
    steady=False
    time=tzero
    
    f=optdens_functions[0](tzero,[0,0,0],0)

    # open file and write dimensions
    file_out=open(filename_optdens, 'w')
    timecell.write2file(file_out,time,True,optdens_tria)

    #cycle in time
    while (time <= tmax) and (not steady ):
        # evalute 
        optdens_tria, steady = timecell.build(
            optdens_tria, steady, 
            optdens_functions,time,bar_cell,flags)
        # write 2 file
        timecell.write2file(file_out,time,False,optdens_tria,steady)
        
        # next time
        time=common.next_time(time)
    file_out.close()
except :
    print( 'No optimal transport density defined ')
    



#######################################################################
# save to vtk
coord_list=[]
topol_list=[]
for i in range(len(coord)):
    coord_list.append([coord[i][0],coord[i][1],0.0])
for i in range(len(topol)):
    topol_list.append(topol[i][0:nnodeincell])


# save all Vtk data for unstructured grid into a vtk object
if (len(optdens_functions) != 0):
    celldata = CellData( Scalars(forcing_tria[:,0],name='forcing'),
                         Scalars(kappa_tria[:,0],name='kappa'),
                         Scalars(tdens0_tria[:,0],name='tdens0'),
                         Scalars(optdens_tria[:,0],name='optdens'))
else :
    celldata = CellData( Scalars(forcing_tria[:,0],name='forcing'),
                         Scalars(kappa_tria[:,0],name='kappa'),
                         Scalars(tdens0_tria[:,0],name='tdens0'))

nodedata = PointData( Scalars(dirac_forcing_nodes[:,0],name='dirac_forcing'))
vtk_object = VtkData( 
    UnstructuredGrid(coord_list,triangle=topol_list),
    celldata,nodedata)
vtk_object.tofile(str(filename_inputs_vtk),'ascii')


##########################
# pflux
scalars.build_and_write(flag_pflux,extra_pflux,
                    tzero,tmax,
                    filename_pflux)

# pmass
scalars.build_and_write(flag_pmass,extra_pmass,
                    tzero,tmax,
                    filename_pmass)

# decay
scalars.build_and_write(flag_decay,extra_decay,
                    tzero,tmax,
                    filename_decay)





#
# build subgrid and subgrid quantities
#
#
# subgrid quanties
#
common.replace('directory','subgrid.fnames',0,str(dir_inputs))
cmd='./subgrid_preprocess/code/subgrid_var.out '
subprocess.call(cmd,shell=True)




