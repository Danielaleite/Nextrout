import networkx as nx
import numpy as np
import itertools
import matplotlib.pyplot as plt
import pickle as pkl
import os
import time
import decimal
#---------------------------------------
import source_sink_generator
import utils
import pre_extraction
#---------------------------------------

def preprocessing_cont(
	folder_name, min_, funct
):
	'''
	Prepocessing the data obtained by the DMK solver to be used as inputs in the pre-extraction step.
	:param folder_name: folder path where the outputs will be stored. It should be written "./runs/folder_name".
	:param min_: threshold for the weights of the edges after filtering.
	:param funct: weights to be assigned to the edges. Either 'tdens' or 'flux'.
	:return:
		G_bar: a weighted networkX graph whose nodes are the barycenters of the elements of the grid. No edges.
		G_triang: the grid graph.
		dict_seq: dictionary mapping grid elements into their vertices. Namely, dict_seq for the key-th element of the
		grid with vertices are n1,n2 and n3 needs to be defined as dict_seq[key]=[n1,n2,n3].  It works with triangular
		and squared grids.
		max_: maximum weight of nodes in G_bar.
	'''

	# Importing the graph structure. We take only the bareycenter and then we compute the edges.

	print('get_baryc')
	start_time = time.time()
	graph_coordinates,_ = utils.get_baryc(folder_name)
	print("--- %s seconds ---" % (time.time() - start_time))

	# Saving the positions of the barycenters in a dict(to plot!)

	print('bar2dict')
	start_time = time.time()
	bar_pos = utils.bar2dict(graph_coordinates)
	print("--- %s seconds ---" % (time.time() - start_time))

	# Extracting the weights (opt_tdens, opt_pot, opt_flux)

	print('extracting_weights')
	start_time = time.time()


	if 'tdens' in funct :
		print(funct)
		opt_tdens = utils.extracting_weights(folder_name,'output/result/'+funct) #<----------------new!
	elif funct == 'flux':
		opt_tdens = utils.extracting_weights(folder_name,'output/result/opt_tdens.dat')
		opt_pot = utils.extracting_weights(folder_name,"output/result/opt_nrm_grad_avg.dat")
	else:
		print('wrong parameter funct:',funct)

	print("--- %s seconds ---" % (time.time() - start_time))

	# Saving the tdens (or pot or flux) in a dict and storing the max

	print('weight2dict')
	start_time = time.time()
	if 'tdens' in funct:
		dict_weights_tdens,weights_tdens = utils.weight2dict(opt_tdens,'output/result/opt_tdens.dat')
	if funct=='flux':
		dict_weights_tdens,weights_tdens= utils.weight2dict(opt_tdens,'output/result/opt_tdens.dat')
		dict_weights_pot,weights_pot= utils.weight2dict(opt_pot,"output/result/opt_nrm_grad_avg.dat")
		dict_weights,weights=utils.prod_dict(dict_weights_tdens,dict_weights_pot,weights_tdens,weights_pot)
	else:
		dict_weights=dict_weights_tdens
		weights=weights_tdens
	max_=max(weights)
	print('The max '+ funct.split('/')[-1].split('.')[0].split('_')[-1]+' in this case is equal to',max_,'. The chosen threshold percentage is', min_*100,'%')
	print('Numerical threshold:',min_*max_)
	print("--- %s seconds ---" % (time.time() - start_time))

	# Completing with zeros

	print('completing_with_zeros')
	start_time = time.time()
	dict_weights = utils.completing_with_zeros(dict_weights,bar_pos)
	print("--- %s seconds ---" % (time.time() - start_time))

	# Defining the graph with the attributes (pos and tdens)

	print('dict2graph')
	start_time = time.time()
	G_bar,pos = pre_extraction.dict2graph(bar_pos,dict_weights)
	print("--- %s seconds ---" % (time.time() - start_time))

	# Getting the grid

	print('getting_grid')
	start_time = time.time()

	file3_ = open(folder_name+"/input/grid.dat", "r")

	grid_ = file3_.readlines()
	file3_.close()
	n_nodes=int(grid_[0])
	grid_nodes=grid_[2:2+n_nodes]
	grid_triangles=grid_[2+n_nodes:]
	#grid_
	print("--- %s seconds ---" % (time.time() - start_time))

	# Saving the grid in a dict

	print("saving_the_grid_to_dict")
	start_time = time.time()
	dict_triang = {}
	n = 1
	pos_triang = [line[:-1].split(" ") for line in grid_nodes]
	for line_ in pos_triang:
		dict_triang[str(n)] = []
		line_ = list(dict.fromkeys(line_))
		if len(line_) < 2:
			line_.append(line_[0])
		for i in line_:
			dict_triang[str(n)].append(float(i))
		n += 1
	for key_ in dict_triang.keys():
		dict_triang[key_] = np.array(dict_triang[key_])
	print("--- %s seconds ---" % (time.time() - start_time))

	# Defining the graph with the nodes of the grid

	print("defining_graph_with_the_nodes_of_the_grid")
	start_time = time.time()
	G_triang = nx.Graph()
	G_triang.add_nodes_from(dict_triang.keys())
	for n, p in dict_triang.items():
		G_triang.nodes[n]["pos"] = p
		# X.nodes[n]['weight'] = weights[int(n)-1]
	pos = nx.get_node_attributes(G_triang, "pos")
	print("--- %s seconds ---" % (time.time() - start_time))

	# Getting the edges of the grid

	print("getting_the_edges_of_the_grid")
	start_time = time.time()
	dict_seq = {}
	n = 1
	pos_ = [line[:-1].split(" ") for line in grid_triangles]
	# print(pos_)
	for line_ in pos_:
		# print(line_)
		dict_seq[str(n)] = []
		line_ = list(dict.fromkeys(line_))
		if len(line_) > 3:
			line_.pop(3)
		for i in line_:
			dict_seq[str(n)].append(str(i))
		n += 1
	for key_ in dict_seq.keys():
		dict_seq[key_] = np.array(dict_seq[key_])

	for i in dict_seq.keys():
		# print(list(itertools.combinations(list(dict_seq[i]), 2)))
		G_triang.add_edges_from(itertools.combinations(list(dict_seq[i]), 2))
	print("--- %s seconds ---" % (time.time() - start_time))

	return G_bar, G_triang, dict_seq
	### end of preprocess


def graph_extraction_from_dat_files(
	folder_path, min_, graph_type, funct, weighting_method,source_flag,sink_flag
):
	'''
	From DMK outputs to pre-extraction step.
	:param folder_path: folder path where the outputs will be stored. It should be written "./runs/folder_name".
	:param min_: threshold for the weights of the edges after pre-extraction.
	:param graph_type: 1 (to use edges and vertices of the grid), 2 (to use only edges), 3 (to use elements of the grid).
	:param funct: weights to be assigned to the edges. Either 'tdens' or 'flux'.
	:param weighting_method: 'ER', 'AVG'.
	:param source_flag: flag used to define the source region in the continuous problem.
    :param sink_flag: flag used to define the sink region in the continuous problem.
	:return:
		Graph: pre-extracted graph.
	'''
	funct_without_dot = funct.split('.')[0]
	new_dir = '../otp_utilities/muffe_sparse_optimization/simplifications/'+folder_path[2:]
	try:
		os.mkdir(new_dir)
	except OSError:
		print ("Creation of the directory %s failed." % new_dir)

	new_dir = '../otp_utilities/muffe_sparse_optimization/simplifications/'+folder_path[2:]+'/'+funct_without_dot
	try:
		os.mkdir(new_dir)
	except OSError:
		print ("Creation of the directory %s failed." % new_dir)


	G_bar, G_triang, dict_seq = preprocessing_cont(
		folder_path+'/', min_, funct
		)

	Graph = pre_extraction.pre_extraction(
		G_bar, G_triang, dict_seq, min_,graph_type,weighting_method
		)

	t_string= '%.0E' % decimal.Decimal( str(min_) )

	print('last part of the path',os.getcwd().split('/')[-1],'=?',folder_path[2:])

	path_='../otp_utilities/muffe_sparse_optimization/simplifications/'+folder_path[2:]+'/'+funct_without_dot+'/'

	file_name = funct_without_dot+'_t'+t_string+'_graph'+str(graph_type)+'_wm'+str(weighting_method)+'.dat'
	print(os.getcwd())
	print('bfore exporting',len(list(nx.connected_components(Graph))))
	with open(path_+file_name, 'wb') as file:
		pkl.dump(Graph, file)


	# Plotting the graph

	print("Plotting.")
	print(os.getcwd())
	###################################################First plot: the extracted graph
	fig, ax = plt.subplots(1, 1, figsize=(25, 25))
	# Frame used to get a plot of the whole space
	frame = nx.path_graph(5)
	pos_frame = {
		0: (0.01, 0.01),
		1: (0.01, 0.99),
		2: (0.99, 0.99),
		3: (0.99, 0.01),
		4: (0.01, 0.01),
	}
	nx.draw(
		frame,
		pos_frame,
		node_size=1,
		width=1,
		edge_color="white",
		node_color="white",
		ax=ax,
	)
	# Plotting the sources and sinks
	if source_flag in [
		"1_rect",
		"3_rect",
		"1_obstacles",
		"3_obstacles",
		"2_rect",
		"2_obstacles",
		"center",
	]:
		ax = source_sink_generator.source_plot(source_flag, ax)
	else:
		source_sink_generator.source_plot(source_flag)
	if sink_flag in ["1_rect", "circle"]:
		ax = source_sink_generator.sink_plot(sink_flag, ax)
	else:
		source_sink_generator.sink_plot(sink_flag)
	# Plotting the simplification
	pos = nx.get_node_attributes(Graph, "pos")
	nx.draw_networkx(
		Graph,
		pos,
		node_size=100,
		with_labels=False,
		edge_color="blue",
		node_color="blue",
		ax=ax,
	)
	# Saving the plots
	t_string = "%.0E" % decimal.Decimal(str(min_))
	file_name = funct_without_dot+'_t'+t_string+'_graph'+str(graph_type)+'_wm'+str(weighting_method)+'.png'

	plt.savefig(path_+file_name)
	plt.show(block=False)
	time.sleep(0.01)
	plt.close("all")

	###################################################Second plot: the extracted graph + colormap based on the btns centrality
	fig, ax = plt.subplots(1, 1, figsize=(25, 25))
	btns_c = nx.betweenness_centrality(Graph)
	nodes = []
	colors = []
	for node in Graph.nodes():
		nodes.append(node)
		colors.append(btns_c[node])

	nx.draw(
		frame,
		pos_frame,
		node_size=1,
		width=1,
		edge_color="white",
		node_color="white",
		ax=ax,
	)
	ec = nx.draw_networkx_edges(Graph, pos, alpha=0.2, ax=ax)
	nc = nx.draw_networkx_nodes(
		Graph,
		pos,
		nodelist=nodes,
		node_color=colors,
		with_labels=False,
		node_size=200,
		cmap=plt.cm.jet,
		ax=ax,
	)
	nc.set_edgecolor("k")

	if source_flag in [
		"1_rect",
		"3_rect",
		"1_obstacles",
		"3_obstacles",
		"2_rect",
		"2_obstacles",
		"center",
	]:
		ax = source_sink_generator.source_plot(source_flag, ax)
	else:
		source_sink_generator.source_plot(source_flag)
	if sink_flag in ["1_rect", "circle"]:
		ax = source_sink_generator.sink_plot(sink_flag, ax)
	else:
		source_sink_generator.sink_plot(sink_flag)

	plt.colorbar(nc)
	plt.axis("on")
	file_name='btns_centr_'+funct_without_dot+'_t'+t_string+'_graph'+str(graph_type)+'_wm'+str(weighting_method)+'.png'
	plt.savefig(path_+file_name)
	plt.show(block=False)
	time.sleep(0.01)
	plt.close("all")

	###################################################Third plot: extracted_graph + colormap given by the opt_funct (tdens or flux)
	fig, ax = plt.subplots(1, 1, figsize=(25, 25))
	nx.draw(
		frame,
		pos_frame,
		node_size=1,
		width=1,
		edge_color="white",
		node_color="white",
		ax=ax,
	)
	edges_, weights_ = zip(
		*nx.get_edge_attributes(Graph, "weight").items()
	)

	ec2 = nx.draw_networkx_edges(
		Graph,
		pos,
		edge_color=weights_,
		width=2.0,
		edge_cmap=plt.cm.jet,
		ax=ax,
	)
	if source_flag in [
		"1_rect",
		"3_rect",
		"1_obstacles",
		"3_obstacles",
		"2_rect",
		"2_obstacles",
		"center",
	]:
		ax = source_sink_generator.source_plot(source_flag, ax)
	else:
		source_sink_generator.source_plot(source_flag)
	if sink_flag in ["1_rect", "circle"]:
		ax = source_sink_generator.sink_plot(sink_flag, ax)
	else:
		source_sink_generator.sink_plot(sink_flag)

	plt.colorbar(ec2)
	plt.axis("on")
	file_name = funct_without_dot+'_distr_t'+t_string+'_graph'+str(graph_type)+'_wm'+str(weighting_method)+'.png'
	plt.savefig(path_+file_name)
	plt.show(block=False)
	time.sleep(0.01)
	plt.close("all")

	#if 'runs/' + os.getcwd().split('/')[-1] == folder_path[2:]:
	#	os.chdir('./../../')

	return Graph


