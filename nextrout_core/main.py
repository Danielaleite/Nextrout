#### import the needed stuff
import dmk_cont
import pre_extraction
import filtering

import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.tri as mtri
import numpy as np
import pickle as pkl



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

    # run the dmk_cont
    grid, subgrid, points, vertices, coord,topol,element_attributes = dmk_cont.grid_gen(ndiv)
    forcing = dmk_cont.forcing_generator(forcing_flag, grid, coord, topol, extra_info=extra_info)
    tdpot = dmk_cont.dmk_cont(forcing,beta_c, ndiv)

    if storing is not None:
        triang = mtri.Triangulation(coord.transpose()[0,:], coord.transpose()[1,:], topol)
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.set_aspect('equal')
        tpc = ax1.tripcolor(triang, -tdpot.tdens,  cmap='gray')
        ax1.tricontour(triang, forcing, levels=40, linewidths=0.1, colors='k')
        fig1.colorbar(tpc)
        plt.savefig(storing+'/dmk_sol.png')
        plt.close()

    # run the graph extraction

    Gpe = pre_extraction.pre_extr(coord, topol, tdpot, min_= min_pe, graph_type=graph_type,weighting_method = weighting_method, DMKw = DMKw)

    if storing is not None:
        weights = np.array([Gpe.edges[edge][DMKw] for edge in Gpe.edges()])
        max_w = max(weights)
        weights/=max_w
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.tricontour(triang, forcing, levels=40, linewidths=0.1, colors='k')
        pos = nx.get_node_attributes(Gpe,'pos')
        nx.draw(Gpe,pos, node_size = 10, node_color = 'k',width = weights*3, ax = ax1)
        plt.savefig(storing+'/Gpe.png')
        plt.close()

    # source/sink generations
    
    sources,sinks = filtering.terminals_from_cont(Gpe, forcing_flag, extra_info, 
        btns_factor_source=btns_factor_source, btns_factor_sink=btns_factor_sink, 
        terminal_criterion=terminal_criterion)

    # run the dmk_discrete

    # get the connected components
    cc_list = list(nx.connected_components(Gpe))

    # apply dmk in the cc
    Gf = nx.Graph()
    count=0
    for cc in cc_list:
        count+=1
        temp_Gpe = Gpe.subgraph(cc)
        temp_sources = [node for node in sources if node in cc]
        temp_sinks = [node for node in sinks if node in cc]

        if len(temp_sources) ==0 or len(temp_sinks) == 0:
            raise ValueError('Not enough sources or sinks.')

        temp_Gf,weights,colors = filtering.filtering(
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
            ax1.tricontour(triang, forcing, levels=40, linewidths=0.1, colors='k')
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
        weights = np.array([Gf.edges[edge][BPw] for edge in Gf.edges()])
        max_w = max(weights)
        weights/=max_w
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.tricontour(triang, forcing, levels=40, linewidths=0.1, colors='k')
        pos = nx.get_node_attributes(Gf,'pos')
        nx.draw(Gf,pos, node_size = 30, node_color = color, width = abs(weights)*3, ax = ax1 )
        plt.savefig(storing+'/Gf.png')
        plt.close()

    
    # storing the results

    if storing is not None:
        if verbose:print('storing at:'+storing)
        files = [('Gpe',Gpe),('Gf',Gf)]
        for ff in files:
            with open(storing + '/'+ff[0]+'.pkl', 'wb') as file:
                pkl.dump(ff[1], file)

    #print('isolated nodes?:',len(list(nx.isolates(Gf))))

    return Gf


### testing:

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

elif forcing_flag == 'rect_cnst2':

    x_source1, y_source1 = (.1,.1)
    wo1 = .2
    ho1 = .5
    rectangles_source = [[(x_source1,y_source1),wo1,ho1]] # bottom left cornner, width, height

    x_sink, y_sink = (.7,.1)
    wi = .2
    hi = .5
    rectangles_sink = [[(x_sink,y_sink),wi,hi]]

    extra_info = [rectangles_source,rectangles_sink]

elif forcing_flag == 'rect_cnst_d':

    x_source1, y_source1 = (.1,.1)
    wo1 = .05
    ho1 = .05
    rectangles_source = [[(x_source1,y_source1),wo1,ho1]] # bottom left cornner, width, height

    x_sink, y_sink = (.7,.1)
    wi = .05
    hi = .05    
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

elif forcing_flag == 'dirac2':

    fplus=[1,2,3,1]
    fminus=[4,2,1]

    xplus=[[0.1,0.2],[0.3,0.4],[0.1,0.7], [0.1,0.9]]
    xminus=[[0.6,0.2],[0.8,0.4],[0.9,0.5]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}

elif forcing_flag == 'dirac3':

    fplus=[3,2]
    fminus=[4,2]

    xplus=[[0.1,0.21],[0.3,0.9]]
    xminus=[[0.6,0.2],[0.8,0.4]]

    Nplus = len(xplus)
    Nminus = len(xminus)

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}

beta_c = 1.5
beta_d = 1.5
flags = ['whole_convex_hull+btns_centr','branch_convex_hull+btns_centr','btns_centr','single']

### running nextrout

nextrout(forcing_flag,
    extra_info,
    beta_c,
    beta_d = beta_d, 
    ndiv = 30, 
    graph_type='1',
    weighting_method = 'ER',
    min_pe = 0.01,
    min_f = 0.1,
    BPw = 'flux',
    weighting_method_simplification ='ER',
    stop_thresh_f = 1e-8,
    verbose = False,
    weight_flag = 'length',
    btns_factor_source=.5, 
    btns_factor_sink=.5,
    terminal_criterion =  flags[3],
    storing = './outputs/')
