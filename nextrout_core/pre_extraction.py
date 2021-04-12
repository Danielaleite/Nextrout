### import needed stuff
import networkx as nx
import numpy as np
from itertools import combinations


def dmk2pre_extr_inputs(coord, topol, tdpot):
	coordinates = {}
	print('len coord 1', len(coord))
	k=-1
	for index_ in coord:
	  k+=1
	  coordinates [k] = index_[:2]

	G_bar = nx.Graph()
	G_triang = nx.Graph()
	k=-1
	centers= {}
	for T in topol: #T is a triangle
	  k+=1
	  edges_in_T = list(combinations(T, 2))
	  for edge in edges_in_T:
	    G_triang.add_edge(edge[0],edge[1])
	  coord_of_nodes_in_T = np.array([coordinates[node] for node in T])
	  center = sum(coord_of_nodes_in_T)/3
	  centers[k] = center
	  G_bar.add_node(k,pos = center, weight = tdpot.tdens[k])

	for node in G_triang.nodes():
	  G_triang.nodes[node]['pos'] = coordinates[node]

	dict_seq = {}
	for key in G_bar.nodes():
	  dict_seq[key] = topol[key]


	return G_bar,G_triang,dict_seq

def get_first_neig(node, dict_seq):
	'''
	This returns the vertices that belong to the triangle for which 'node' is the barycenter
	:param node: node in G_bar.
	:param dict_seq:  dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:return:
		same_triang: all the nodes in the triangle (or element of the grid) s.t., node is its barycenter.
	'''
	same_triang = list(dict_seq[node])
	return same_triang



def get_sec_neig(node, dict_seq):
	'''
	This returns all the triangles that share either a vertex or an edge with the triangle in which the 'node'
	is the barycenter.
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(
		node, dict_seq
	)  # getting the nodes of the triang that has 'node' as barycenter
	dict_sec_neig = (
		{}
	)  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []

	# indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles but also their bar
	index_ = {}
	index_[node] = []
	for key in dict_seq.keys():
		# if one of the nodes is in the triangle, then ...
		if (
				same_triang[0] in dict_seq[key]
				or same_triang[1] in dict_seq[key]
				or same_triang[2] in dict_seq[key]
		):

			# (to avoid repeating the triangles)
			repeated_tr = list(np.array(same_triang) == np.array(dict_seq[key]))

			if repeated_tr != [True, True, True]:
				# ... we define the triangle indexes for a particular node
				dict_sec_neig[node].append(dict_seq[key])
				index_[node].append(key)
	return dict_sec_neig, index_



def get_sec_neig_edges(node, dict_seq):
	'''
	This returns all the triangles that share an edge with the triangle in which the 'node' is the barycenter
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(
		node, dict_seq
	)  # getting the nodes of the triangle that has 'node' as barycenter
	dict_sec_neig = (
		{}
	)  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []
	index_ = (
		{}
	)  # indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles
	# but also their bar
	index_[node] = []
	for key in dict_seq.keys():
		if (
				(same_triang[0] in dict_seq[key] and same_triang[1] in dict_seq[key])
				or (same_triang[1] in dict_seq[key] and same_triang[2] in dict_seq[key])
				or (same_triang[2] in dict_seq[key] and same_triang[0] in dict_seq[key])
		):
			repeated_tr = list(np.array(same_triang) == np.array(dict_seq[key]))
			if repeated_tr != [True, True, True]:
				dict_sec_neig[node].append(dict_seq[key])
				index_[node].append(key)
	return dict_sec_neig, index_

def get_sec_neig_edges_square(node, dict_seq, node2box_index):
	'''
	This returns all the triangles that share an edge with the triangle in which the 'node' is the barycenter
	:param node: target node.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
	:return:
		dict_sec_neig: dictionary, s.t., dict_sec_neig[node]= indexes of all the surrounding triangles of 'node'.
		index_: **check this output. it seems to be removable.**
	'''
	same_triang = get_first_neig(node, dict_seq)  # getting the nodes of the triang that has 'node' as barycenter
	neighboring_boxes = []
	pair_of_vertices_same_triang = list(itertools.combinations(same_triang, 2))
	for pair in pair_of_vertices_same_triang:
		node1=pair[0]
		node2=pair[1]
		indexes_node1 = node2box_index[node1]
		indexes_node2 = node2box_index[node2]
		print(indexes_node1,indexes_node2)
		intersect = [index for index in indexes_node1 if index in indexes_node2]
		neighboring_boxes = neighboring_boxes + intersect
	neighboring_boxes = list(set(neighboring_boxes))
	neighboring_boxes.remove(node)
	# print(neighboring_boxes)
	dict_sec_neig = {}  # dict that contains all the surrounding triangles for the 'node'
	dict_sec_neig[node] = []

	# indexes of all the surrounding triangles. This will be useful to call not only the mentioned triangles but also their bar
	index_ = neighboring_boxes

	return dict_sec_neig, index_







def connecting_edges(G_bar, node, min_, graph_type, dict_seq, max_, weighting_method, input_flag=None,node2box_index=None):
	'''
	Testing the condition of the tdens for a single barycenter . If condition is satisfied, add the edge with weight
	equal to the min. In-place modification of G_bar.
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	:param node: target node.
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges), 3 (to use elements of the grid).
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:param max_: maximum weight of nodes in G_bar.
	:param weighting_method: 'ER', 'AVG'.
	:param input_flag: 'image' or None (for dat files)
	:param node2box_index:  given a node, node2box_index[node] returns the indices of all the elements of the partition s.t.,
	node is a vertex of these elements
	:return:
	'''
	if graph_type == "1":
		# If True, then just the 'first neighbors' are taken.
		if input_flag != 'image':
			#node = str(node)
			dict_sec_neig, index_ = get_sec_neig(node, dict_seq)
			index_ = index_[node]
		else:
			node=node-1
			dict_sec_neig, index_ = get_sec_neig_square(node, dict_seq, node2box_index)

	elif graph_type == "2":
		# If True, then just the 'second neighbors' are taken (firsts included).
		if input_flag != 'image':
			node = node
			dict_sec_neig, index_ = get_sec_neig_edges(node, dict_seq)
			index_ = index_[node]
		else:
			node = node - 1
			dict_sec_neig, index_ = get_sec_neig_edges_square(node, dict_seq, node2box_index)


	#print(index_)
	for bar_ in index_:
		if graph_type == "1" or graph_type == "2":  # both edges and vertices
			# Checking if the weights are greater than the threshold x max
			if (
					G_bar.nodes[bar_]["weight"] > min_ * max_
					and G_bar.nodes[node]["weight"] > min_ * max_
			):
				if weighting_method == "AVG":
					# If True, then we assign weights to the edges based on the 'AVG' method
					w_list = [G_bar.nodes[bar_]["weight"], G_bar.nodes[node]["weight"]]
					w = sum(w_list) / len(w_list)  # <--- the average
					# w=w0 #<--- this is the minimum
					G_bar.add_edge(node, bar_, weight=w)
				elif weighting_method == "ER":
					# If True, then we just add the edge. Later on we define the weight.
					G_bar.add_edge(node, bar_)


def grid_filter(G_bar, G_triang, min_, dict_seq):
	'''
    This script adds the edges of a triangle to the graph type 3 if the weight of its barycenter is greater than the
    threshold (min_) x max_.
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	:param G_triang: the grid graph.
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:return:
		 G_pre_extracted: the input graph but with edges generated according to graph_type 3.
	'''

	nodes_dict = G_bar.nodes(data='weight')
	max_ = max([entry[1] for entry in nodes_dict])

	G_pre_extracted = G_triang.copy()

	# Getting the max key (why?)

	#keys = [int(key) for key in G_pre_extracted.nodes().keys()]
	#max_label = max(keys)

	# Removing the edges of the triangulation graph. We are left only with the nodes. Also the useless ones.

	edges_triang = list(G_pre_extracted.edges())
	G_pre_extracted.remove_edges_from(edges_triang)

	for bar_ in G_bar.nodes():

		# Getting the edges of the triangle for a given bar 'bar_'

		same_triang = get_first_neig(bar_, dict_seq)
		w = G_bar.nodes[bar_]["weight"]

		if w >= min_ * max_:
			# If True, we add the edges to the graph
			for node in same_triang:
				G_pre_extracted.nodes[node]['weight']=0
			edges = list(itertools.combinations(same_triang, 2))
			membership = [
				edges[0] in G_pre_extracted.edges(),
				edges[1] in G_pre_extracted.edges(),
				edges[2] in G_pre_extracted.edges(),
			]

			# We iterate over them to define their weights
			for i in [0, 1, 2]:
				bool_val = membership[i]
				# If the edge is in the graph (thanks to another barycenter), we sum the weight of the barycenter 'bar_' to the existing weight
				if bool_val == True:
					new_weight = (
							G_pre_extracted.edges[(edges[i][0], edges[i][1])]["weight"]
							+ float(w) / 2.0
					)
					G_pre_extracted.edges[(edges[i][0], edges[i][1])]["weight"] = new_weight
				# If the edge is not in the graph, we add half of the weight of the barycenter
				else:
					new_weight = float(w) / 2.0
					G_pre_extracted.add_edge(edges[i][0], edges[i][1], weight=new_weight)

	return G_pre_extracted


def node_edge_filter(G_bar, min_, graph_type, dict_seq, weighting_method,input_flag=None,node2box_index=None):
	'''
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
	:param min_:  threshold for the weights of the edges after pre-extraction.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges).
	:param dict_seq:  dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:param weighting_method: 'ER', 'AVG'.
	:param input_flag: 'image' or None (for dat files)
	:param node2box_index: given a node, node2box_index[node] returns the indices of all the elements of the partition s.t.,
	node is a vertex of these elements
	:return:
		G_bar: the input graph but with edges generated according to graph_type.
	'''
	"""
	This script generates graph type 1 and 2.
	"""
	nodes_dict = G_bar.nodes(data='weight')
	max_ = max([entry[1] for entry in nodes_dict])

	reduced_node_list = [int(node) for node in G_bar.nodes() if G_bar.nodes[node]['weight']>min_ * max_]

	# Iterate over all the numbers (<--> nodes) to test the condition about the threshold:
	for n in reduced_node_list:#range(len(G_bar.nodes()))
		connecting_edges(
			G_bar, n, min_, graph_type, dict_seq, max_, weighting_method,input_flag,node2box_index
		)

	if weighting_method == "ER":

		# Compute degree centrality

		deg = nx.degree_centrality(G_bar)
		N = len(G_bar.nodes())

		# If True, then apply 'ER' method

		for edge in G_bar.edges():
			G_bar.edges[(edge[0], edge[1])]["weight"] = G_bar.nodes[edge[0]]["weight"] / (
					deg[edge[0]] * (N - 1)
			) + G_bar.nodes[edge[1]]["weight"] / (deg[edge[1]] * (N - 1))
	elif weighting_method == "AVG":
		pass
	else:
		print("weighting_method not defined")

	return G_bar


"""
Filtering triangles
"""


def graph_builder(
		G_bar, G_triang, dict_seq, min_, graph_type='1',weighting_method='ER'
):
	'''
	:param G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid.
	:param G_triang: the grid graph.
	:param dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular and squared grids.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges), 3 (to use elements of the grid).
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param max_: maximum weight of nodes in G_bar.
	:param weighting_method: 'ER', 'AVG'.
	:return:
		G_pre_extracted: pre-extracted graph.
	'''


	if graph_type == "3" and weighting_method == "ER":
		print("this wm does not apply for graph type 3. Try with 'AVG'.")
	else:

		if graph_type == "1" or graph_type == "2":
			G_pre_extracted = G_bar.copy()
			# filtering
			edges_ = list(G_pre_extracted.edges())
			G_pre_extracted.remove_edges_from(edges_)
			G_pre_extracted = node_edge_filter(
				G_pre_extracted, min_, graph_type, dict_seq, weighting_method
			)
		elif graph_type == "3":
			# filtering
			# print('here')
			G_pre_extracted = grid_filter(G_bar, G_triang, min_, dict_seq)

		"""
		Removing the non-useful barycenters
		"""

		G_pre_extracted.remove_nodes_from(list(nx.isolates(G_pre_extracted)))

		return G_pre_extracted

def pre_extr(coord, topol, tdpot, min_):

	G_bar,G_triang,dict_seq = dmk2pre_extr_inputs(coord, topol, tdpot)

	Gpe = graph_builder(G_bar,G_triang,dict_seq, min_) 

	return Gpe