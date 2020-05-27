import time
import networkx as nx
import numpy as np
import operator
import itertools
#--------------------
import quality_measure
#--------------------

def pixel_circle_counter(center,R, partition_dict, G):
    intersecting_circle=[]
    for key in G.nodes():
        corners = partition_dict[key+1]
        bool_list = [(center[0]-R<=point[0]<=center[0]+R) & (center[1]-R<=point[1]<=center[1]+R) for point in corners]
        if 0<sum(bool_list)<4:
            intersecting_circle.append(key+1)
    return intersecting_circle


def single_box_finder(intersection_list, partition_dict):  # under construction
    # placing the first square in the first box
    box = [intersection_list[0]]
    box_coord = partition_dict[intersection_list[0]]
    # print('bc',box_coord)
    for key in intersection_list[1:]:
        # testing if there's an intersection
        # print(key)
        chosen_coord = partition_dict[key]
        # print('cc',chosen_coord)
        A = len(box_coord)
        # print('a',A)
        # print('--u',chosen_coord+box_coord)
        # print('--su',set(chosen_coord+box_coord))
        AunionB = len(set(chosen_coord + box_coord))
        # print('ab',AunionB)
        if AunionB != A + 4:
            box_coord = list(set(box_coord + chosen_coord))
            box = box + [key]
    return box, box_coord

def cluster_counter(intersection_list,partition_dict):
    box={}
    len_intersect = len(intersection_list)
    box_total=0
    while len_intersect>0:
        box,_=single_box_finder(intersection_list,partition_dict)
        intersection_list = [n for n in intersection_list if n not in box]
        len_intersect=len(intersection_list)
        box_total+=1
    return box_total

def branch_counter(intersection_list,partition_dict):
    box_total = cluster_counter(intersection_list,partition_dict)
    if box_total==2 or box_total==0:
        #path
        nbr = 0
    elif box_total==1:
        #endpoint
        nbr=-1
    else:
        #bifurcation,trifurcations,...
        nbr=1
    return nbr


def terminal_finder(D, partition_dict, G):
    print('graph size:' ,len(G.nodes) ,' nodes' ,len(G.edges) ,' edges.')
    start_time = time.time()
    terminal_list =[]
    color_nbr =[]
    nodes_for_correction = {}
    nbr_graph =nx.Graph()

    # betweenness centrality
    bn =nx.betweenness_centrality(G)

    # time
    print("bn centralities: --- %s second(s) ---" % (time.time() - start_time))
    start_time = time.time()

    # closeness centrality
    cn =nx.closeness_centrality(G)

    # time
    print("cn centralities: --- %s second(s) ---" % (time.time() - start_time))
    start_time = time.time()

    # building the centers of the masks
    number_of_masks =round( 1 /(D))
    print('number of masks:' , 1 /D ,number_of_masks)
    x_coord_mask =list(np.linspace( D /2 , 1 - D /2 ,number_of_masks))
    x_y_coord_mask =list(itertools.product(x_coord_mask ,x_coord_mask))
    # print(x_y_coord_mask)
    # sweeping the masks
    N=0
    for center in x_y_coord_mask:
        N+=1
        # print('N',N)
        intersection_list = pixel_circle_counter(center,D/2, partition_dict, G)
        nbr = branch_counter(intersection_list,partition_dict)

        nodes_in_the_mask=[node for node in G.nodes()
                           if (
                             center[0]-D/2<= G .nodes[node]['pos'][0]<center[0]+D/2)
                           & (center[1]-D/2<= G .nodes[node]['pos'][1]<center[1]+D/2)
                           ]
        nbr_graph.add_nodes_from(nodes_in_the_mask)

        subGraph=G.subgraph(nodes_in_the_mask)

        if nbr==-1 and len(nodes_in_the_mask)>0:
            # the least betweenness centrality
            bn_in_the_mask = {n:cn[n] for n in nodes_in_the_mask}
            term = min(bn_in_the_mask.items(), key=operator.itemgetter(1))[0]
            terminal_list.append(term)
        elif nbr>0 and len(nodes_in_the_mask)>0:
            # get the connected components and then
            cn_in_the_mask = {n:bn[n] for n in nodes_in_the_mask}
            term = max(cn_in_the_mask.items(), key=operator.itemgetter(1))[0]
            # terminal_list.append(term)

            # compute the closeness center for each one

        # add this node to the terminal l is t

        elif nbr==0:
            ccom = list(nx.connected_components(subGraph))
            lcc = len(ccom )

            if lcc==1:
                # print('nothing special in this 0-branch')
                pass
            elif lcc>1:
                nodes_for_correction[N]=[ ]

                Dn=D*(1.1)
                nodes_in_the_extended_mask=[node for node in G.nodes()
                                            if (center[0]-Dn/2<=G. nodes[node]['pos'][0]<center[0]+Dn/2)
                                            & (center[1]-Dn/2<=G. nodes[node]['pos'][1]<center[1]+Dn/2)
                                            ]


                subGraph_extended =G.subgraph(nodes_in_the_extended_mask)
                ccom_extended = list(nx.connected_components(subGraph_extended))
                lcc_e = len(ccom_extended)
                if True:  # lcc==lcc_e:

                    for elem in ccom:
                        elem_subgraph = subGraph.subgraph(elem)
                        bn_sub = {n:bn[n] for n in elem_subgraph.nodes()}
                        term = min(bn_sub.items(), key=operator.itemgetter(1))[0]
                        terminal_list.append(term)

                        # correction
                        nodes_for_correction[N].append(term)


        for node in nodes_in_the_mask:
            nbr_graph.nodes[node]['pos'] = G.nodes[node]['pos']
            if nbr == -1:
                color_nbr.append('red')
            elif nbr == 0:
                color_nbr.append('blue')
            else:
                color_nbr.append('green')

    print("rest: --- %s second(s) ---" % (time.time() - start_time))
    return nbr_graph, color_nbr, terminal_list, nodes_for_correction, number_of_masks