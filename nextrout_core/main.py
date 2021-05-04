#### import the needed stuff
import dmk_cont
import pre_extraction
import filtering

import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.tri as mtri
import numpy as np
import pickle as pkl
import os



def nextrout(
    forcing_flag,
    extra_info,
    beta_c,
    beta_d = 1,
    ndiv=15, 
    graph_type='1',
    weighting_method = 'ER',
    DMKw = 'tdens',
    min_pe = 0.01, 
    min_f = 0.001,
    verbose = True, 
    weighting_method_simplification = 'ER',
    BPw = 'tdens0', 
    stop_thresh_f = 1e-6, 
    btns_factor_source=.1, 
    btns_factor_sink=.1,
    terminal_criterion = 'btns_centr',
    weight_flag = 'unit',
    storing = None
    ):
    
    inputs = {}
    inputs['continuous'] = {}
    # run the dmk_cont
    grid, subgrid, points, vertices, coord,topol,element_attributes = dmk_cont.grid_gen(ndiv)
    forcing, triang_source_indices,triang_sink_indices = dmk_cont.forcing_generator(forcing_flag, grid, coord, topol, extra_info=extra_info)
    tdpot, timefun = dmk_cont.dmk_cont(forcing,beta_c, ndiv, storing = storing)

    print(os.getcwd())

    if storing is not None:
        triang = mtri.Triangulation(coord.transpose()[0,:], coord.transpose()[1,:], topol)
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.set_aspect('equal')
        tpc = ax1.tripcolor(triang, -tdpot.tdens,  cmap='gray')
        ax1.tricontour(triang, forcing, cmap='RdBu_r')
        fig1.colorbar(tpc)
        plt.savefig(storing+'/dmk_sol.png')
        plt.close()

        inputs['continuous']['ndiv'] = ndiv
        inputs['continuous']['forcing'] = forcing
        inputs['continuous']['pflux'] = beta_c

    # run the graph extraction

    tdens_weights = tdpot.tdens

    Gpe = pre_extraction.pre_extr(coord, topol, tdens_weights,triang_source_indices,triang_sink_indices, min_= min_pe, graph_type=graph_type,weighting_method = weighting_method, DMKw = DMKw)

    node_colors = []
    for node in Gpe.nodes():
    	terminal_val = Gpe.nodes[node]['terminal']
    	if terminal_val == 1:
    		node_colors.append('g')
    	elif terminal_val == -1:
    		node_colors.append('r')
    	else:
    		node_colors.append('k')

    if storing is not None:
        weights = np.array([Gpe.edges[edge][DMKw] for edge in Gpe.edges()])
        max_w = max(weights)
        weights/=max_w
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.tricontour(triang, forcing, cmap='RdBu_r')
        pos_Gpe = nx.get_node_attributes(Gpe,'pos')
        nx.draw(Gpe,pos_Gpe, node_size = 10, node_color = node_colors, width = weights*3, ax = ax1)
        plt.savefig(storing+'/Gpe.png')
        plt.close()

    # source/sink generations
    
    sources,sinks = filtering.terminals_from_cont(Gpe, forcing_flag, extra_info, 
        btns_factor_source=btns_factor_source, btns_factor_sink=btns_factor_sink, 
        terminal_criterion=terminal_criterion)

    # run the dmk_discrete
    print('len',len(sources),len(sinks))
    # get the connected components
    cc_list = list(nx.connected_components(Gpe))
    print('CC',len(cc_list))
    # apply dmk in the cc
    Gf = nx.Graph()
    count=0
    for cc in cc_list:
        count+=1
        temp_Gpe = Gpe.subgraph(cc)
        temp_sources = [node for node in sources if node in cc]
        temp_sinks = [node for node in sinks if node in cc]

        if len(temp_sources) ==0 or len(temp_sinks) == 0:
            raise ValueError('Not enough sources or sinks. Increase btns_factor.')

        temp_Gf,weights,colors, inputs_discr = filtering.filtering(
            temp_Gpe, 
            temp_sources, 
            temp_sinks, 
            beta_d = beta_d, 
            tdens0 = 2, # 2 means not unitary (i.e., taken from Gpe)
            threshold = min_f, 
            BPweights = BPw, 
            stopping_threshold_f = stop_thresh_f)
        
        # put everything together

        Gf = nx.disjoint_union(Gf, temp_Gf)

        if storing is not None and len(cc_list)!=1:
            fig1, ax1 = plt.subplots(figsize=(10, 10))
            ax1.tricontour(triang, forcing, cmap='RdBu_r')
            pos = nx.get_node_attributes(temp_Gf,'pos')
            nx.draw(temp_Gf,pos, node_size = 30, node_color = colors, width = abs(weights)*2 ,ax = ax1)
            plt.savefig(storing+'/Gf_'+str(count)+'.png')
            plt.close()

    

    deg = nx.degree_centrality(Gf)

    if weighting_method_simplification == 'ER':

        N = len(Gf.nodes())

        for edge in Gf.edges():

            Gf.edges[(edge[0], edge[1])][BPw] = Gf.nodes[edge[0]][
                                                                 'weight'] / (
                                                                     deg[edge[0]] * (N - 1)) + \
                                                             Gf.nodes[edge[1]][
                                                                 'weight'] / (
                                                                     deg[edge[1]] * (N - 1))
    elif weighting_method_simplification == 'IBP':
        
        raise ValueError('not implemented yet.')
        '''
        Gf_relabeled = relabeling(Gf, Gpe)
        for edge in Gf_relabeled.edges():
            Gf.edges[edge[0],edge[1]][BPw]=Gpe.edges[edge[0],edge[1]]['weight']
        '''
        
    elif weighting_method_simplification == 'BPW':
        pass

    

    if storing is not None:
        if len(cc_list)==1:
            color = colors
        else:
            color = 'k'
        weights = np.array([abs(Gf.edges[edge][BPw]) for edge in Gf.edges()])
        max_w = max(weights)
        weights/=max_w
        print('max',max_w)

        edge_labels = {}
        for edge in Gf.edges():
            edge_labels[edge]=round(abs(Gf.edges[edge][BPw])/max_w,2)
            if edge_labels[edge] == 0:
                edge_labels[edge] = 0

        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.tricontour(triang, forcing, cmap='RdBu_r')
        pos = nx.get_node_attributes(Gf,'pos')
        nx.draw(Gf,pos, node_size = 30, node_color = color, width = abs(weights)*3, ax = ax1 )
        nx.draw_networkx_edge_labels(Gf,pos,edge_labels=edge_labels,font_color='red', font_size = 8, ax = ax1)
        plt.savefig(storing+'/Gf.png')
        plt.close()

    
    # storing the results

    if storing is not None:
        if verbose:print('storing at:'+storing)
        files = [('Gpe',Gpe),('Gf',Gf)]
        for ff in files:
            with open(storing + '/'+ff[0]+'.pkl', 'wb') as file:
                pkl.dump(ff[1], file)

        # storing inputs
        inputs['discrete'] = inputs_discr

        with open(storing + '/inputs.pkl', 'wb') as file:
                pkl.dump(inputs, file)

        pos = nx.get_node_attributes(Gf,'pos')

        with open(storing + '/Gf_node_locations.pkl', 'wb') as file:
                pkl.dump(pos, file)

    #print('isolated nodes?:',len(list(nx.isolates(Gf))))

    return Gf


