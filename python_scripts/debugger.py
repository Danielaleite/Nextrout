import os
import pickle as pkl
import shutil
import distutils
import numpy as np
import shutil
import distutils
from threading import Thread
import re
import networkx as nx
import matplotlib.pyplot as plt
#------------------------------------
import continuous2graph
import discrete2graph
import filtering
import pre_extraction
import quality_measure
import utils
import source_sink_generator
#------------------------------------

network_extraction_script = False
extraction_from_image_script = True
quality_measure_script = True

what_to_test = [network_extraction_script, extraction_from_image_script, quality_measure_script]

if network_extraction_script:

    data_flags = ['10_b12_10dv_sf15_rect_cnstrect_cnst']#,
                  #'parabola_1_b10_20dv_sf15_5rch5rch' ]#'delta_1_b12_18dv_sf15_rect_cnstrect_cnst']#,                           ]
    cont_data={}
    disc_data={}
    differences={}
    G_={}
    G_['new']={}
    G_['known']={}

    flags = {data_flags[0]:[
        '10',
    '1.2',
    '10',
    'rect_cnst',
    'rect_cnst',
    '1.5',
    'yes',
    'yes',
    'yes,no'
    ]}

    '''
    ,
        data_flags[1]: [
        'parabola_1',
    '1.0',
    '20',
    '5rch',
    '5rch',
    '1.5',
    'yes',
    'yes',
    'yes,no'
    ]
    }
    '''


    for folder in data_flags:

        cont_data[folder]={}
        disc_data[folder]={}
        differences[folder]={}

        flag_list = flags[folder][0]
        beta_list = flags[folder][1]
        ndiv_list = flags[folder][2]
        source_list = flags[folder][3]
        sink_list = flags[folder][4]
        bd_list = flags[folder][5]
        dmk_input = flags[folder][6]
        ge_input = flags[folder][7]
        gs_input = flags[folder][8]

        '''
        Read the well know information
        '''

        # Read the number of files in the folder

        cont_data[folder]['known']= [
            os.path.join('./debugger_input_files/Continuous/'+folder,file)
            for file in os.listdir('./debugger_input_files/Continuous/'+folder)
        ]

        for thing in cont_data[folder]['known']:
            #isFile = os.path.isfile(thing)
            if thing.split('/')[-1]=='input':
                nbr_input_files_known = 0
                l = len([
                    os.path.join(thing,file)
                    for file in os.listdir(thing)])
                nbr_input_files_known+=l
                print('files in the input folder:',nbr_input_files_known)
            elif thing.split('/')[-1]=='output':
                nbr_output_files_known = 0
                l = len([
                    os.path.join(thing+'/result', file)
                    for file in os.listdir(thing+'/result')])
                nbr_output_files_known += l
                print('files in the output/result folder:', nbr_output_files_known)

        disc_data[folder]['known']= [
            os.path.join('./debugger_input_files/Discrete/'+folder+'/tdens',file)
            for file in os.listdir('./debugger_input_files/Discrete/'+folder+'/tdens')
        ]


        nbr_discrete_files_known= len(disc_data[folder]['known'])
        print('files in the discrete folder:',nbr_discrete_files_known)

        #Load the graphs
        known_graphs_in_disc_data = [
            os.path.join('./debugger_input_files/Discrete/' + folder + '/tdens', file)
            for file in os.listdir('./debugger_input_files/Discrete/' + folder + '/tdens') if '.dat' in file
        ]

        for graph in known_graphs_in_disc_data:
            if len(graph) > 30:  # taking  a look at the simpl
                with open(graph, 'rb') as file:
                    G_['known'][folder]= pkl.load(file)

        '''
          Run a new simulation
          '''



        if sink_list == '=':
            sink_list = source_list
        source_list = source_list.split(',')
        sink_list = sink_list.split(',')
        source_sink_list = [[source_list[i], sink_list[i]] for i in range(len(source_list))]
        print(source_sink_list)
        for flag in flag_list.split(','):
            for beta in beta_list.split(','):
                for ndiv in ndiv_list.split(','):
                    for source_sink in source_sink_list:
                        for beta_discr in bd_list.split(','):
                            folder_name = "%s" % flag + "_" + "b" + str(
                                int(10 * (float(beta)))) + '_' + "%s" % ndiv + "dv" + '_sf' + str(
                                int(10 * (float(beta_discr)))) + '_' + "%s" % source_sink[0] + "%s" % source_sink[1]

                            if dmk_input == 'yes':
                                #click.echo('Your .ctrl file is ready! ')
                                # click.echo('Your flag is %s!' %flag)
                                file = open('inputs.ctrl', 'w+')

                                file.write("# DOMAIN ####################################################### " + "\n")
                                file.write('2d ! flag_domain (2d, 3d, surface)' + "\n")

                                file.write('# MESH #########################################################' + "\n")
                                file.write('rect_cnst path/grid.dat ! flag_grid' + "\n")
                                file.write('%s' % ndiv + ' ! ndiv' + "\n")
                                file.write('1 ' + ' ! nref' + "\n")  # --------------------------- !!!!!!!!!!

                                file.write('# PREBUILD EXAMPLES ############################################' + "\n")
                                file.write('%s' % beta + ' extra_path ! flag_pflux ' + "\n")
                                file.write('1.0 extra_path ! flag_pmass (gamma)' + "\n")
                                file.write('1.0 extra_path ! flag_decay' + "\n")
                                file.write('1.0 extra ! decay0' + "\n")
                                file.write('rect_cnst  path/frog_source.dat ! flag_source' + "\n")  # %s' %source_sink[0]+
                                file.write('rect_cnst  path/frog_sink2.dat ! flag_sink' + "\n")  # %s' %source_sink[1]+
                                file.write('0 ! flag_normalize' + "\n")
                                # file.write('%s' %flag + '!flag_tdens0 ' + "\n" )
                                file.write('%s' % flag + ' ' + 'path/tdens0.dat ! flag_tdens0 ' + "\n")
                                file.write('1.0 extra_path ! flag_kappa' + "\n")

                                file.write('## TIME RANGE ###################################################' + "\n")
                                file.write('0.0 ! tzero' + "\n")
                                file.write('5.0e2 ! tmax (stop while loop)')
                                # list = file.readlines()
                                file.close()


                                # folder_name =  "%s" %flag + "_" + "b"  +str(int(10*(float(beta)))) + '_' + "%s" %ndiv + "dv" + '_sf'+str(int(10*(float(beta_discr))))+ '_'+ "%s" %source_sink[0]+ "%s" %source_sink[1]

                                def init(folder_name):
                                    new_dir = "./runs/" + folder_name
                                    try:
                                        os.mkdir(new_dir)
                                    except OSError:
                                        print("Creation of the directory %s failed." % new_dir)


                                try:
                                    shutil.rmtree("./runs/" + folder_name)
                                except OSError:
                                    pass
                                command = "./dmk_folder.py assembly " + './runs/' + folder_name + " inputs.ctrl "
                                os.system(command)

                                source_sink_generator.source_sink_generator('./runs/' + folder_name, ndiv, source_sink[0], source_sink[1])
                                if source_sink[0] != 'rect_cnst' and source_sink[1] != 'rect_cnst':
                                    source_sink_preprocess('./runs/' + folder_name)

                                command = "./dmk_folder.py run " + './runs/' + folder_name + " new_muffa.ctrl > outputs_dmk_c.txt"
                                os.system(command)

                                command = "./dmk_folder.py get-graph " + './runs/' + folder_name + " 0.1 " + " 100000000000  > outputs_gg.txt"
                                os.system(command)

                                command = "./dmk_folder.py vtk -a -tdens " + './runs/' + folder_name + " > outputs_vtk.txt"
                                os.system(command)


                            else:
                                print('Skipping DMK-solver part.')
                            ############## GRAPH EXTRACTION######################

                            errors = []
                            for funct in ['tdens']:
                                print('Now we are running: ', funct)

                                for graph_type in ['1']:  # ,'2','3']:
                                    print('Now we are running graph: ', graph_type)
                                    for weighting_method_graph in ['ER']:  # , 'AVG']:
                                        # print('Now we are running weighting method graph: ', weighting_method_graph)
                                        if funct == 'tdens':
                                            if 0 < float(beta) <= 1:
                                                t_list = [.1]  # ,.65]
                                            else:
                                                t_list = [0.001]  # <------------ !!
                                        else:  # flux
                                            if 0 < float(beta) <= 1:
                                                t_list = [.1]  # ,.65]
                                            else:
                                                t_list = [0.001]

                                        for threshold in t_list:
                                            subfolder = './runs/' + folder_name
                                            t = float(threshold)

                                            new_dir = subfolder + '/' + funct

                                            try:
                                                os.mkdir(new_dir)
                                            except OSError:
                                                print("Creation of the directory %s failed." % new_dir)

                                            if ge_input == 'yes':
                                                print('Flag and tolerance:', subfolder, t, graph_type, funct)
                                                G = continuous2graph.graph_extraction_from_dat_files(subfolder, t, graph_type,funct, weighting_method_graph,source_sink[0],source_sink[1])
                                            else:
                                                print('Skipping graph-extraction part.')
                                            print(gs_input.split(','))
                                            if gs_input.split(',')[0] == 'yes':
                                                for minimum in [0.001]:  # [0.001]:
                                                    for weighting_method_simplification in ['ER']:  # IBP, BPW
                                                        for btns_factor in [(1, 1)]:  # 0.010]
                                                            print('=======================================================',
                                                                  graph_type, funct, weighting_method_graph,
                                                                  weighting_method_simplification)
                                                            min_ = float(minimum)
                                                            btns_factor_source = float(btns_factor[0])
                                                            btns_factor_sink = float(btns_factor[1])
                                                            BP_weights = 'BPtdens'
                                                            reduction_flag = gs_input.split(',')[
                                                                1]  # write 'yes' to get the 2nd simpl
                                                            # print('>>>>>____________________Computing bp simplification for',subfolder, funct, t, graph_type, min_,weighting_method_graph, weighting_method_simplification)
                                                            i = 0
                                                            #beta_discr = float(beta_discr)
                                                            i += 1
                                                            #updating_beta_discrete(beta_discr)
                                                            # try:
                                                            Simp, conv_report = discrete2graph.graph_filtering_from_dat_files(subfolder, t, graph_type,
                                                                                                  beta_discr,
                                                                                                  funct, min_,
                                                                                     btns_factor_source, btns_factor_sink,
                                                                                     weighting_method_graph,
                                                                                     weighting_method_simplification,
                                                                                     source_sink[0], source_sink[1],
                                                                                     BP_weights, reduction_flag)



        #os.chdir('../../Tests/')
        # Read the number of files in the folder

        '''
        Read the well know information
        '''

        # Read the number of files in the folder

        cont_data[folder]['new'] = [
            os.path.join('./debugger_input_files/Continuous/' + folder, file)
            for file in os.listdir('./debugger_input_files/Continuous/' + folder)
        ]

        for thing in cont_data[folder]['new']:
            # isFile = os.path.isfile(thing)
            if thing.split('/')[-1] == 'input':
                nbr_input_files_new = 0
                l = len([
                    os.path.join(thing, file)
                    for file in os.listdir(thing)])
                nbr_input_files_new += l
                print('files in the input folder:', nbr_input_files_new)
            elif thing.split('/')[-1] == 'output':
                nbr_output_files_new = 0
                l = len([
                    os.path.join(thing + '/result', file)
                    for file in os.listdir(thing + '/result')])
                nbr_output_files_new += l
                print('files in the output/result folder:', nbr_output_files_new)

        disc_data[folder]['new'] = [
            os.path.join('./debugger_input_files/Discrete/' + folder + '/tdens', file)
            for file in os.listdir('./debugger_input_files/Discrete/' + folder + '/tdens')
        ]

        nbr_discrete_files_new = len(disc_data[folder]['new'])
        print('files in the discrete folder:', nbr_discrete_files_new)

        # Load the graphs
        known_graphs_in_disc_data = [
            os.path.join('./debugger_input_files/Discrete/' + folder + '/tdens', file)
            for file in os.listdir('./debugger_input_files/Discrete/' + folder + '/tdens') if '.dat' in file
        ]

        differences[folder]=[nbr_input_files_known-nbr_input_files_new,nbr_output_files_known-nbr_output_files_new,nbr_discrete_files_known-nbr_discrete_files_new]

        for graph in known_graphs_in_disc_data:
            if len(graph)>30:#taking  a look at the simpl
                with open(graph, 'rb') as file:
                    G_['new'][folder] = pkl.load(file)

        '''
        Compare the new results with the know ones
        '''
        #print(G_)

        # print number of nodes, number of edges of both cases

        #Compute the difference between the vectors generated by both graphs on a specific attribute
    print('_______________________________________________')
    print('Simulations finished! Time to test the results:')
    for folder in data_flags:
        print('\nTested folder:',folder)
        print('Testing number of files in folders:')
        print('Difference in the number of FILES between the known data and the new one (cont/input/,cont/output/result/,discr/):',differences[folder])
        g_known=G_['known'][folder]
        #print(g_known.nodes(data=True))
        g_new=G_['new'][folder]
        #print(g_new.nodes(data=True))
        print('Difference in the number of NODES in known graph and the new one:',len(g_known.nodes)-len(g_new.nodes))
        print('Difference in the number of EDGES in known graph and the new one:',len(g_known.edges)-len(g_new.edges))
        #node attributes
        pos=np.array([g_known.nodes[node]['pos'] for node in g_known.nodes()])-np.array([g_new.nodes[node]['pos'] for node in g_new.nodes()])
        print('Sum of the difference in the POSITIONS (node-wise) in known graph and the new one:',sum(pos) )
        w_node = np.array([g_known.nodes[node]['weight'] for node in g_known.nodes()]) - np.array([g_new.nodes[node]['weight'] for node in g_new.nodes()])
        print('Sum of the difference in the WEIGHTS (node-wise) in known graph and the new one:', sum(w_node))
        pot = np.array([g_known.nodes[node]['op_pot'] for node in g_known.nodes()]) - np.array([g_new.nodes[node]['op_pot'] for node in g_new.nodes()])
        print('Sum of the difference in the OPT_POT (node-wise) in known graph and the new one:', sum(pot))
        terminal = np.array([g_known.nodes[node]['terminal'] for node in g_known.nodes()]) - np.array([g_new.nodes[node]['terminal'] for node in g_new.nodes()])
        print('Sum of the difference in the TERMINAL PROPERTY (S+,S-) (node-wise) in known graph and the new one:', sum(terminal))
        #edge attributes
        w_edge=np.array([g_known.edges[edge]['weight'] for edge in g_known.edges()])-np.array([g_new.edges[edge]['weight'] for edge in g_new.edges()])
        print('Sum of the difference in the WEIGHTS (edge-wise) in known graph and the new one:',sum(w_edge) )
        print('Convergence_report:',conv_report)


if extraction_from_image_script:
    ### Testing Image Extraction ###
    G_image={}
    G_image['known']={}
    G_image['new']={}

    new_size=100
    ratio=new_size/1200
    #print('ratio:',ratio)
    t1 = 0.1
    t2 = .5
    image_path = "./runs/graph_from_image/image.jpg"
    number_of_cc=1
    beta_d = 1
    number_of_colors=50
    graph_type='1'

    path_to_graph =  "./debugger_input_files/Images/image/"


    for graph in ['extracted_graph.pkl','filtered_graph.pkl']:
        path = path_to_graph + graph
        with open(path, 'rb') as file:
            G_image['known'][graph]= pkl.load(file)

    '''
    Run a new simulation
    '''

    filtering.img2filtering(image_path, new_size, number_of_colors, t1, t2, number_of_cc, graph_type, beta_d)

    path_to_graph =  "./runs/image/"

    for graph in ['extracted_graph.pkl','filtered_graph.pkl']:
        path = path_to_graph + graph
        with open(path, 'rb') as file:
            G_image['new'][graph]= pkl.load(file)


    '''
    Comparing
    '''

    for graph in ['extracted_graph.pkl','filtered_graph.pkl']:
        print('\nTested folder:', path_to_graph)
        print('Tested graph:',graph)
        g_known=G_image['known'][graph]
        #print(g_known.nodes(data=True))
        g_new=G_image['new'][graph]
        #print(g_new.nodes(data=True))
        node_diff = len(g_known.nodes) - len(g_new.nodes)
        edge_diff = len(g_known.edges)-len(g_new.edges)
        print('Difference in the number of NODES in known graph and the new one:',node_diff)
        print('Difference in the number of EDGES in known graph and the new one:',edge_diff)
        #node attributes
        if node_diff == 0 :
            pos=np.array([g_known.nodes[node]['pos'] for node in g_known.nodes()])-np.array([g_new.nodes[node]['pos'] for node in g_new.nodes()])
            print('Sum of the difference in the POSITIONS (node-wise) in known graph and the new one:',sum(pos) )
            w_node = np.array([g_known.nodes[node]['weight'] for node in g_known.nodes()]) - np.array([g_new.nodes[node]['weight'] for node in g_new.nodes()])
            print('Sum of the difference in the WEIGHTS (node-wise) in known graph and the new one:', sum(w_node))
            if 'extr' not in graph:
                terminal = np.array([g_known.nodes[node]['terminal'] for node in g_known.nodes()]) - np.array([g_new.nodes[node]['terminal'] for node in g_new.nodes()])
                print('Sum of the difference in the TERMINAL PROPERTY (S+,S-) (node-wise) in known graph and the new one:', sum(terminal))
        #edge attributes
        elif edge_diff == 0:
            w_edge=np.array([g_known.edges[edge]['weight'] for edge in g_known.edges()])-np.array([g_new.edges[edge]['weight'] for edge in g_new.edges()])
            print('Sum of the difference in the WEIGHTS (edge-wise) in known graph and the new one:',sum(w_edge) )
        else:
            print('edge and node properties not compared.')

if quality_measure_script:
    '''
    Load data
    '''
    path = "./debugger_input_files/Quality_measure/"
    with open(path+'./qw_3-29.pkl', 'rb') as file:
        qm_weights_old = pkl.load(file)
    with open(path+'./L_3-29.pkl', 'rb') as file:
        qm_L_old = pkl.load( file)

    '''
    Run simulation
    '''
    qm = {}
    # print('\n      N    | Pre-extracted Graph    |')
    # print('           |       qw   |  L        |')
    for N in range(3, 30):
        qm[N] = {}

        partition_dict, dict_seq, node2box_index = quality_measure.partition_set(N)
        _, G_triang = quality_measure.partition(N)

        colors = utils.horizontal_line(N)

        G_bar = pre_extraction.weighted_partition2bar_graph(partition_dict, colors)

        G_pre_extracted = G_bar.copy()
        # print(G_bar.nodes(data=True))
        # filtering
        edges_ = list(G_pre_extracted.edges())
        G_pre_extracted.remove_edges_from(edges_)
        # print('getting graph')

        graph_type = "1"
        min_ = .5
        metric = "l1"

        G_pre_extracted = pre_extraction.node_edge_filter(G_pre_extracted,
                                                          .5,
                                                          graph_type,
                                                          dict_seq,
                                                          'ER',
                                                          'image',
                                                          node2box_index)  # 12

        q_measure, q_weight, L, triang_weight_vector, weight_difference_vector = quality_measure.q_measure(G_pre_extracted,
                                                                                                           G_bar,
                                                                                                           partition_dict,
                                                                                                           metric, min_=min_)
        qm[N]['qw'] = q_weight
        qm[N]['L'] = L

        # print("%10.3e" % N, '|', "%f" % qm[N]['qw'], '  |', "%f " % qm[N]['L'], '|')

    qm_weights_new = np.array([qm[N]['qw'] for N in qm.keys()])
    qm_L_new = np.array([qm[N]['L'] for N in qm.keys()])



    '''
    Compare data
    '''
    print('\nTested folder:', path)

    len_qw_diff = len(qm_weights_new) - len(qm_weights_old)
    len_L_diff = len(qm_L_old) - len(qm_L_old)
    print('Size difference for qw:', len_qw_diff)
    print('Size difference for L:', len_L_diff)
    if len_qw_diff==0 and len_L_diff==0:
        diff_qw = sum(qm_weights_new - qm_weights_old)
        diff_L = sum(qm_L_new - qm_L_old)
        print('Substracting the arrays (weight):',diff_qw)
        print('Substracting the arrays (L):', diff_L)
    else:
        print('arrays do not have the same size.')
