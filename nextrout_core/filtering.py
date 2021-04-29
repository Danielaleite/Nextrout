## import stuff
import networkx as nx
import numpy as np
import sys

from scipy.spatial import ConvexHull, convex_hull_plot_2d
from scipy.spatial import distance

#root = '/home/dleite/community_detection_nextrout/Nextrout'
root = '/home/dtheuerkauf/Nextrout'

# Import I/O for timedata
try:
    sys.path.append(root+'/dmk_utilities/globals/python/timedata/')
    import timedata as td
except:
    print("Global repo non found")

# Import geometry tools
sys.path.append(root+'/dmk_utilities/geometry/python/')
import meshtools as mt
sys.path.append(root+'/dmk_utilities/dmk_solver/otp_solver/preprocess/assembly/')
import example_grid

# Import dmk tools
sys.path.append(root+'/dmk_utilities/dmk_solver/otp_solver/python/')
import dmk_p1p0 
sys.path.append(root+'/dmk_utilities/dmk_solver/build/python/fortran_python_interface/')
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
sys.path.append(root+'/dmk_utilities/dmk_solver/graph_otp_solver/python')
import dmk_graph




def concatenate(lists):
    '''
    This concatenates all the lists contained in a list.
    :param lists: a list of lists.
    :return:
        new_list: concatenated list.
    '''
    new_list = []
    for i in lists:
        new_list.extend(i)
    return new_list


def terminals_from_cont(Graph, forcing_flag, extra_info, btns_factor_source, btns_factor_sink, terminal_criterion='branch_convex_hull+btns_centr'):
    '''
    Computation of source and sink nodes. This script uses information about the inputs of the DMK solver.
    There are three criteria for the selection of the nodes.
    :param Graph: a networkx graph to be filtered.
    :param source_flag: flag used to define the source region in the continuous problem.
    :param sink_flag: flag used to define the sink region in the continuous problem.
    :param btns_factor_source: threshold for the nodes in the source region (see more in paper).
    :param btns_factor_sink: threshold for the nodes in the sink region (see more in paper).
    :param terminal_criterion: 'branch_convex_hull+btns_centr' (combination of btns centr and convex hull by branches),
    'whole_convex_hull+btns_centr' (combination of btns centr and convex hull of the source and sink regions),
    'btns_centr' (only btns centr).

    :return:
        possible_terminals_source: for each i, possible_terminals_source[i]= "sources" of i-th cc.
        possible_terminals_sink: for each i, possible_terminals_sink[i]= "sources" of i-th cc.
    '''

    # Defining source-sink nodes of the graph (nodes inside the source or sink of the continuous DMK)

    nodes_in_source = []
    nodes_in_sink = []

    for node in Graph.nodes():
        terminal_val = Graph.nodes[node]['terminal']
        if terminal_val == 1:
            nodes_in_source.append(node)
        elif terminal_val == -1:
            nodes_in_sink.append(node)

    '''
    if 'rect' in forcing_flag:
        
        rectangles_source = extra_info[0] # it should be a list of lists: every sublist has 3 elements: (x,y), w, h
        rectangles_sink = extra_info[1]

        for node in Graph.nodes():

          x,y = Graph.nodes[node]['pos']

          for rect in rectangles_source:

            x_source = rect[0][0]
            y_source = rect[0][1]
            wo = rect[1]
            ho = rect[2]

            if (x_source <= x <= x_source+wo) and (y_source <= y <= y_source+ho):
              nodes_in_source.append(node)

          for rect in rectangles_sink:

            x_sink = rect[0][0]
            y_sink = rect[0][1]
            wi = rect[1]
            hi = rect[2]

            if (x_sink <= x <= x_sink+wi) and (y_sink <= y <= y_sink+hi):
              nodes_in_sink.append(node)


    elif 'dirac' in forcing_flag:
        
        Nplus = extra_info['Nplus']
        Nminus = extra_info['Nminus']

        xplus = extra_info['xplus'] 
        xminus = extra_info['xminus']
        
        nodes_in_source = []
        nodes_in_sink = []

        for node in Graph.nodes():
            for xy in xplus:
                xy = np.array(xy)
                if distance.euclidean(Graph.nodes[node]['pos'],xy)<.1:
                    nodes_in_source.append(node)
            for xy in xminus:
                xy = np.array(xy)
                if distance.euclidean(Graph.nodes[node]['pos'],xy)<.1:
                    nodes_in_sink.append(node)
    else:
        raise ValueError('forcing flag not recognized.')
    '''
    bn = nx.betweenness_centrality(Graph, normalized=True)

    # min bn inside the source and sink

    max_bn_source = max([bn[node] for node in nodes_in_source])
    max_bn_sink = max([bn[node] for node in nodes_in_sink])

    # Defining source-sink candidates based only on btns

    kind_of_leaf_nodes_source = [key for key in nodes_in_source if
                                 bn[key] <= btns_factor_source * max_bn_source]  # *(min_bn_source)+.0001]
    kind_of_leaf_nodes_sink = [key for key in nodes_in_sink if
                               bn[key] <= btns_factor_sink * max_bn_sink]  # *(min_bn_sink+.0001)]

    # Defining the subgraphs induced by the candidates

    sub_source = Graph.subgraph(kind_of_leaf_nodes_source)
    sub_sink = Graph.subgraph(kind_of_leaf_nodes_sink)

    # Removing repeated labels

    possible_terminals_source = set(kind_of_leaf_nodes_source)
    possible_terminals_sink = set(kind_of_leaf_nodes_sink)

    print('poss',possible_terminals_source,possible_terminals_sink)
    
    if terminal_criterion == 'single':
        
        terminals = list(possible_terminals_source)+list(possible_terminals_sink)
        
        #possible_terminals_source = list(possible_terminals_source)
        
        extra_nodes = []
        
        for i in range(len(terminals)-1):#range(len(possible_terminals_source)-1):
            
            
            node1 = terminals[i]
            
            pos1 = Graph.nodes[node1]['pos']
            
            stop = False
            
            
            for j in range(i+1,len(terminals)):
                
                node2 = terminals[j]
                
                if node1 != node2 and not stop:
                    
                    pos2 = Graph.nodes[node2]['pos']
                    dst = distance.euclidean(pos1, pos2)
                    
                    if dst < .1:
                        
                        # check distance
                        
                        if node1 not in extra_nodes:
                            
                            stop = True
                        
                            extra_nodes.append(node1)
            
        possible_terminals_source = [int(node) for node in possible_terminals_source if node not in extra_nodes]
        possible_terminals_sink = [int(node) for node in possible_terminals_sink if node not in extra_nodes]


    elif terminal_criterion != 'btns_centr':
        # Getting the coordinates to compute convex hulls

        coordinates_in_source = np.asarray(
            [[Graph.nodes[node]['pos'][0], Graph.nodes[node]['pos'][1]] for node in
             nodes_in_source])
        coordinates_in_source_list = concatenate(list(coordinates_in_source))

        coordinates_in_sink = np.asarray(
            [[Graph.nodes[node]['pos'][0], Graph.nodes[node]['pos'][1]] for node in
             nodes_in_sink])
        coordinates_in_sink_list = concatenate(list(coordinates_in_sink))

        # If the number of coordinates (or nodes) is not more than 3, then no convex hull computation

        if len(coordinates_in_source) >= 4 and len(coordinates_in_sink) >= 4:

            # Computing convex hull for the branches

            if terminal_criterion == 'branch_convex_hull+btns_centr':
                source_hull = ConvexHull(coordinates_in_source)
                index_source_hull = np.asarray(source_hull.vertices)
                nodes_source_hull = []
                coord_source_hull = np.asarray(
                    [coordinates_in_source[node] for node in index_source_hull])
                for node in nodes_in_source:
                    x, y = Graph.nodes[node]['pos']
                    if x in coord_source_hull[:, 0] and y in coord_source_hull[:, 1]:
                        nodes_source_hull.append(node)

                sink_hull = ConvexHull(coordinates_in_sink)
                index_sink_hull = np.asarray(sink_hull.vertices)
                nodes_sink_hull = []
                coord_sink_hull = np.asarray([coordinates_in_sink[node] for node in index_sink_hull])
                for node in nodes_in_sink:
                    x, y = Graph.nodes[node]['pos']
                    if x in coord_sink_hull[:, 0] and y in coord_sink_hull[:, 1]:
                        nodes_sink_hull.append(node)
                possible_terminals_source = set(kind_of_leaf_nodes_source + nodes_source_hull)
                possible_terminals_sink = set(kind_of_leaf_nodes_sink + nodes_sink_hull)

            # Computing convex hull for all the nodes defined as candidates

            elif terminal_criterion == 'whole_convex_hull+btns_centr':  # not working!
                single_source_hull = ConvexHull(coordinates_in_source_list)
                single_index_source_hull = np.asarray(single_source_hull.vertices)
                nodes_source_hull = []
                single_coord_source_hull = np.asarray(
                    [coordinates_in_source_list[node] for node in single_index_source_hull])
                for node in nodes_in_source:
                    x, y = Graph.nodes[node]['pos']
                    if x in single_coord_source_hull[:, 0] and y in single_coord_source_hull[:, 1]:
                        nodes_source_hull.append(node)
                single_sink_hull = ConvexHull(coordinates_in_sink_list)
                single_index_sink_hull = np.asarray(single_sink_hull.vertices)
                nodes_sink_hull = []
                single_coord_sink_hull = np.asarray(
                    [coordinates_in_sink_list[node] for node in single_index_sink_hull])
                for node in nodes_in_sink:
                    x, y = Graph.nodes[node]['pos']
                    if x in single_coord_sink_hull[:, 0] and y in single_coord_sink_hull[:, 1]:
                        nodes_sink_hull.append(node)
                possible_terminals_source = set(kind_of_leaf_nodes_source + nodes_source_hull)
                possible_terminals_sink = set(kind_of_leaf_nodes_sink + nodes_sink_hull)


    return possible_terminals_source, possible_terminals_sink



def filtering(Gpe, sources = None, sinks = None,beta_d = 1.5,threshold=1e-3, tdens0 = 1, BPweights = 'tdens', stopping_threshold_f = 1e-6, weight_flag='unit', rhs = None):


    if sources is None and sinks is None and rhs is None:

        raise ValueError('Either rhs or sources/sinks need to be passed as inputs.')



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

    # tdens0

    if tdens0 != 1:
        try:
            tdens0 = np.array([(Gpe_rel.edges[edge]['tdens']) for edge in edges])
        except:
            tdens0 = np.array([(Gpe_rel.edges[edge]['flux']) for edge in edges])

    # topol

    topol =np.zeros((nedges,2))
    k=-1
    for edge in edges:
      k+=1
      topol[k,:] = edge 

    # weight (uniform)

    weight =np.empty(nedges, dtype=object)
    
    k=-1
    for edge in edges:
        k+=1
        if weight_flag == 'unit':
            weight[k] = 1
        elif weight_flag == 'length':
            weight[k] = distance.euclidean(Gpe.nodes[edge[0]]['pos'],Gpe.nodes[edge[1]]['pos'])

    print('weights',weight)
    # rhs (f+ and f-)

    if sinks is not None and sources is not None: # there are lists from the sources and sinks are going to be chosen.
    # (else) if this is not pass, then the rhs is passed.

        rhs = np.zeros(nnodes)
        sources_rel = [mapping[node] for node in sources]
        sinks_rel = [mapping[node] for node in sinks]

        number_sources = len(sources_rel)
        number_sinks = len(sinks_rel)

        for node in nodes:
            if node in sources_rel:
                rhs[node] = 1/number_sources
            elif node in sinks_rel:
                rhs[node] = -1/number_sinks
            else:
                rhs[node] = 0
    else:
        sources_rel = [i for i in range(len(rhs)) if rhs[i]>0]
        sinks_rel = [i for i in range(len(rhs)) if rhs[i]<0]

    assert sum(rhs) < .01
    assert len(rhs) == nnodes
    print('rhs',rhs)
    # init and set controls
    ctrl = Dmkcontrols.DmkCtrl()
    Dmkcontrols.get_from_file(ctrl,root+'/nextrout_core/dmk_discr.ctrl')
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

    [info,tdens,pot,flux,timefun] = dmk_graph.dmk_graph(
        topol,
        rhs,
        pflux = beta_d,
        #tdens0 =  tdens0,
        tolerance = stopping_threshold_f,
        weight= weight,
        ctrl = ctrl)

    if (info==0):
        print('Convergence achieved')

    max_flux = max(flux)
    max_tdens = max(tdens)
    Gf = nx.Graph()
    ed_count = -1
    weights_in_Gf = []
    for edge in Gpe_rel.edges():
        ed_count+=1
        if BPweights == 'flux':
            if abs(flux[ed_count])> max_flux*threshold:
                Gf.add_edge(*edge,  flux = flux[ed_count])
                weights_in_Gf.append(flux[ed_count])

        elif BPweights == 'tdens':
            if abs(tdens[ed_count])> max_tdens*threshold:
                Gf.add_edge(*edge,  tdens = tdens[ed_count])
                weights_in_Gf.append(tdens[ed_count])

        else:
            raise ValueError('BPweights flag not defined!.')
        try:
            Gf.add_node(edge[0], weight = Gpe_rel.nodes[edge[0]]['tdens'])# todo: this needs to be fixed once the flux is working again (BPweights)
            Gf.add_node(edge[1], weight = Gpe_rel.nodes[edge[1]]['tdens'])
        except:
            pass


    Gf.remove_nodes_from(list(nx.isolates(Gf)))
        
    weights_in_Gf = np.array(weights_in_Gf)
    colors = []

    for node in Gf.nodes():

        Gf.nodes[node]['pos'] = Gpe_rel.nodes[node]['pos']

        if node in sources_rel:
            colors.append('g')
        elif node in sinks_rel:
            colors.append('r')
        else:
            colors.append('k')

    return Gf, weights_in_Gf, colors