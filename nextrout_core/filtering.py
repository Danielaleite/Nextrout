## import stuff
import networkx as nx
import numpy as np
import sys

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
sys.path.append('../dmk_utilities/dmk_solver/graph_otp_solver/python')
import dmk_graph


def filtering(Gpe, threshold=0.0001):

	### relabeling

	mapping = {}
	k=-1
	for node in Gpe.nodes():
	  k+=1
	  mapping[node] = k
	Gpe_rel  =nx.relabel_nodes(Gpe, mapping, copy=True)


	edges = Gpe_rel.edges()
	nedges = len(edges)
	nodes = Gpe_rel.nodes()
	nnodes = len(nodes)
	leaves = [node for node in Gpe_rel.nodes() if Gpe_rel.degree(node)==1]

	# topol

	topol =np.zeros((nedges,2))
	k=-1
	for edge in edges:
	  k+=1
	  topol[k,:] = edge
	print('topol',topol[:5])  

	# weight (uniform)

	weight =np.empty(nedges, dtype=object)
	k=-1
	for edge in edges:
	  k+=1
	  weight[k] = 1
	print('weight',weight[:5])

	# rhs (f+ and f-)
	import random
	rhs = np.zeros(nnodes)
	sources = leaves[:1]
	number_sources = len(sources)
	number_sinks = len(leaves)-len(sources)
	for node in nodes:
	  if node in leaves:
	    if node in sources:
	      rhs[node] = 1/number_sources
	    else:
	      rhs[node] = -1/number_sinks
	  else:
	    rhs[node] = 0
	print('balanced?',sum(rhs))


	# init and set controls
	ctrl = Dmkcontrols.DmkCtrl()
	Dmkcontrols.get_from_file(ctrl,'../dmk_utilities/dmk_solver/graph_otp_solver/python/examples/FaccaBenzi2021_TestCase1/dmk.ctrl')
	# if and where save data
	ctrl.id_save_dat=1
	ctrl.fn_tdens='tdens.dat'
	ctrl.fn_pot='pot.dat'
	# if and where save log
	ctrl.id_save_statistics=1
	ctrl.fn_statistics='dmk.log'
	# if print info 
	# 
	print(ctrl.outer_solver_approach)

	[info,tdens,pot,flux,timefun] = dmk_graph.dmk_graph(topol,rhs,1.0,1e-6,weight,ctrl)
	if (info==0):
	    print('Convergence achieved')


	Gf = nx.Graph()
	ed_count = -1
	for edge in Gpe_rel.edges():
		ed_count+=1
		if flux[ed_count]> threshold:
			Gf.add_edge(*edge,  flux = flux[ed_count])
	print('len of Gf:',len(Gf.nodes()))
	return Gf