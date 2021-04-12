#### import the needed stuff
import dmk_cont
import pre_extraction
import filtering

import matplotlib.pyplot as plt
import networkx as nx
import matplotlib.tri as mtri
import numpy as np




def nextrout(forcing_flag,extra_info,beta_c,beta_d,ndiv=15, min_pe = 0.01, min_f = 0.001,verbose = True):

    # run the dmk_cont
    grid, subgrid, points, vertices, coord,topol,element_attributes = dmk_cont.grid_gen(ndiv)
    forcing = dmk_cont.forcing_generator(forcing_flag, grid, coord, topol, extra_info=extra_info)
    tdpot = dmk_cont.dmk_cont(forcing,beta_c, ndiv)

    if verbose:
        triang = mtri.Triangulation(coord.transpose()[0,:], coord.transpose()[1,:], topol)
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.set_aspect('equal')
        tpc = ax1.tripcolor(triang, tdpot.tdens,  cmap='GnBu')
        #fig1.colorbar(tpc)
        #tpc = ax1.tripcolor(triang, dmkin.penalty_weight,   cmap='cubehelix')
        fig1.colorbar(tpc)
        plt.show()

    # run the graph extraction
    Gpe = pre_extraction.pre_extr(coord, topol, tdpot, min_= min_pe)

    if verbose:
        weights = np.array([Gpe.edges[edge]['weight'] for edge in Gpe.edges()])
        max_w = max(weights)
        weights/=max_w
        fig1, ax1 = plt.subplots(figsize=(10, 10))
        ax1.tricontour(triang, forcing, levels=40, linewidths=0.1, colors='k')
        pos = nx.get_node_attributes(Gpe,'pos')
        nx.draw(Gpe,pos, node_size = 10, node_color = 'k',width = weights,ax = ax1)
        plt.show()
    
    # source/sink generations
    flags = ['whole_convex_hull+btns_centr','branch_convex_hull+btns_centr','btns_centr','single']
    sources,sinks = filtering.terminals_from_cont(Gpe, forcing_flag, extra_info, 
        btns_factor_source=.1, btns_factor_sink=.1, 
        terminal_criterion=flags[3])

    # run the dmk_discrete

    Gf,weights,colors = filtering.filtering(Gpe, sources, sinks, beta_d = beta_d, threshold = min_f, BPweights = 'flux')

    if verbose:

        pos = nx.get_node_attributes(Gf,'pos')
        nx.draw(Gf,pos, node_size = 30, node_color = colors, width = abs(weights)*5 )
        plt.show()

    # storing the results!

    ### pending

    return Gf


### testing:

# generate forcing

forcing_flag = 'dirac'

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

elif forcing_flag == 'dirac2':
    Nplus=2
    Nminus=2

    fplus=[3,2]
    fminus=[4,2]

    xplus=[[0.1,0.21],[0.3,0.9]]
    xminus=[[0.6,0.2],[0.8,0.4]]

    extra_info = {'Nplus':Nplus,
                   'Nminus':Nminus,
                    'fplus':fplus,
                    'fminus':fminus,
                    'xplus':xplus,
                    'xminus':xminus}

beta_c = 1.5
beta_d = 1.5
nextrout(forcing_flag,extra_info,beta_c,beta_d, ndiv = 20, min_f = 0.001)
