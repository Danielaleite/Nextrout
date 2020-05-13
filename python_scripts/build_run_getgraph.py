#!/usr/bin/env python

import click
import os
import shutil
import distutils
import numpy as np
import shutil
import distutils
from threading import Thread
from continuous2graph import *
from discrete2graph import *
from filtering import *
from source_sink_generator import *
import re
import networkx as nx
import pickle as pkl
import matplotlib.pyplot as plt



import os, os.path
import shutil
import glob
import subprocess
from shutil import copytree, ignore_patterns


click.echo('Define the parameters for the continuous part here!')
@click.command()
@click.option('--flag_list', prompt='What is your flag?', help='This is where you should put the name of the flag.')
@click.option('--ndiv_list', prompt='What is the number of divisions for the Mesh?', help='This is the number of divisions for your triangular mesh.')
@click.option('--beta_list', prompt='What is the value for beta?', help='This is the value of your Beta exponent.')
@click.option('--source_list', prompt='What is(are) the flag(s) for source?', help='Source information.')
@click.option('--sink_list', prompt='What is(are) the flag(s) for sink?', help='Sink information.')
@click.option('--bd_list', prompt='What is the val for the discrete beta?', help='Beta discrete information.')
@click.option('--dmk_input', prompt='continuous DMK solver?', help='Execute the DMK solver.')
@click.option('--ge_input', prompt='graph extraction?', help='Extract the graph from the DMK solutions.')
@click.option('--gs_input', prompt='graph simplification, reduction?', help='Simplify the graph using discrete DMK solver. Do a second-step simplification.')

def new(flag_list,beta_list,ndiv_list,source_list,sink_list, bd_list,dmk_input,ge_input,gs_input):
	#beta_discr=1.5#---------------------------------------------------------------------------------------------------!
	if sink_list == '=':
		sink_list=source_list
	source_list = source_list.split(',')
	sink_list = sink_list.split(',')
	source_sink_list = [ [source_list[i],sink_list[i] ] for i in range(len(source_list))]
	print(source_sink_list)
	for flag in flag_list.split(','):
		for beta in beta_list.split(','):
			for ndiv in ndiv_list.split(','):
				for source_sink in source_sink_list:
					for beta_discr in bd_list.split(','):
						folder_name =  "%s" %flag + "_" + "b"  +str(int(10*(float(beta)))) + '_' + "%s" %ndiv + "dv" + '_sf'+str(int(10*(float(beta_discr))))+ '_'+ "%s" %source_sink[0]+ "%s" %source_sink[1]
							
						if dmk_input == 'yes':
							click.echo('Your .ctrl file is ready! ')
							#click.echo('Your flag is %s!' %flag)
							file=open('inputs.ctrl','w+')

							file.write("# DOMAIN ####################################################### " + "\n")
							file.write('2d ! flag_domain (2d, 3d, surface)' + "\n")

							file.write('# MESH #########################################################' + "\n")
							file.write('rect_cnst path/grid.dat ! flag_grid' + "\n")
							file.write('%s' %ndiv + ' ! ndiv' + "\n")
							file.write('0 ' + ' ! nref' + "\n")#--------------------------- !!!!!!!!!!

							file.write('# PREBUILD EXAMPLES ############################################' + "\n")
							file.write('%s' %beta + ' extra_path ! flag_pflux ' + "\n")
							file.write('1.0 extra_path ! flag_pmass (gamma)' + "\n")
							file.write('1.0 extra_path ! flag_decay' + "\n")
							file.write('1.0 extra ! decay0' + "\n")
							file.write('rect_cnst  path/frog_source.dat ! flag_source' + "\n")#%s' %source_sink[0]+
							file.write('rect_cnst  path/frog_sink2.dat ! flag_sink' + "\n")#%s' %source_sink[1]+
							file.write('0 ! flag_normalize' + "\n")
							#file.write('%s' %flag + '!flag_tdens0 ' + "\n" )
							file.write('%s' %flag + ' ' + 'path/tdens0.dat ! flag_tdens0 ' + "\n" )
							file.write('1.0 extra_path ! flag_kappa' + "\n")

							file.write('## TIME RANGE ###################################################' + "\n")
							file.write('0.0 ! tzero' + "\n")
							file.write('5.0e2 ! tmax (stop while loop)')
							#list = file.readlines()
							file.close()
							
							#folder_name =  "%s" %flag + "_" + "b"  +str(int(10*(float(beta)))) + '_' + "%s" %ndiv + "dv" + '_sf'+str(int(10*(float(beta_discr))))+ '_'+ "%s" %source_sink[0]+ "%s" %source_sink[1]
							
							def init(folder_name):
								new_dir="./runs/"+folder_name
								try:
									os.mkdir(new_dir)
								except OSError: 
									print("Creation of the directory %s failed." % new_dir)
							try:
								shutil.rmtree("./runs/"+folder_name)
							except OSError:
								pass	
							command="./dmk_folder.py assembly " + './runs/'+folder_name + " inputs.ctrl "
							os.system(command)

							source_sink_generator('./runs/'+folder_name, ndiv, source_sink[0], source_sink[1])
							if source_sink[0]!='rect_cnst' and source_sink[1]!='rect_cnst':
								source_sink_preprocess('./runs/'+folder_name)

							command="./dmk_folder.py run " + './runs/'+ folder_name  + " new_muffa.ctrl > outputs_dmk_c.txt"
							os.system(command)

							command="./dmk_folder.py get-graph " +  './runs/'+ folder_name  + " 0.1 " + " 100000000000  > outputs_gg.txt"
							os.system(command)

							command="./dmk_folder.py vtk -a -tdens " +  './runs/'+ folder_name +" > outputs_vtk.txt"
							os.system(command)
							
							

							'''
							#moving the rhs.dat so we can the perturbation
							

							print(os.getcwd())
							os.system('cp ./runs/'+folder_name+'/input/rhs.dat ' + '../globals/python_timedata/rhs.dat')
							os.chdir('../globals/python_timedata/')
							os.system('python perturbation.py')
							os.system('python balance.py rhs_perturbation_damped.dat rhs_integrated.dat')
							os.chdir('../../Tests')
							os.system('cp ../globals/python_timedata/rhs_integrated.dat'+'  ./runs/'+folder_name+'/input/rhs_integrated.dat ')
							
							
							
							#running dmk once more
							
							command="./dmk_folder.py run " + './runs/'+ folder_name  + " new_muffa.ctrl > outputs_dmk_c.txt"
							os.system(command)
							command="./dmk_folder.py vtk -a -tdens " +  './runs/'+ folder_name +" > outputs_vtk.txt"
							os.system(command)
							'''
							#shutil.copyfile('inputs.ctrl', 'Test_set/inputs.ctrl')

						else:
							print('Skipping DMK-solver part.')
						############## GRAPH EXTRACTION######################
						
						errors=[]
						for funct in ['tdens']:
							print('Now we are running: ', funct)

							for graph_type in ['1']:#,'2','3']:
								print('Now we are running graph: ', graph_type)
								for weighting_method_graph in ['ER']:#, 'AVG']:
									#print('Now we are running weighting method graph: ', weighting_method_graph)
									if funct =='tdens':
										if 0<float(beta)<=1:
											t_list=[.1]#,.65]
										else:
											t_list=[0.001]#<------------ !!
									else: #flux
										if 0<float(beta)<=1:
											t_list=[.1]#,.65]
										else:
											t_list=[0.001]

									for threshold in t_list:
										subfolder='./runs/'+folder_name
										t=float(threshold)

										new_dir = subfolder+'/'+funct

										try:
											os.mkdir(new_dir)
										except OSError:
											print ("Creation of the directory %s failed." % new_dir)
										
										
										if ge_input == 'yes':
											print('Flag and tolerance:',subfolder,t,graph_type,funct)
											G=graph_extraction_from_dat_files(subfolder, t, graph_type,funct, weighting_method_graph,source_sink[0],source_sink[1]) 
										else:
											print('Skipping graph-extraction part.')
										print(gs_input.split(','))
										if gs_input.split(',')[0]=='yes':
											for minimum in [0.001]:#[0.001]:
												for weighting_method_simplification in ['ER']:#IBP, BPW
													for btns_factor in [(-1,-1)]:# 0.010]
														print('=======================================================',graph_type,funct,weighting_method_graph,weighting_method_simplification)
														min_=float(minimum)
														btns_factor_source=float(btns_factor[0])
														btns_factor_sink=float(btns_factor[1])
														BP_weights='BPtdens'
														reduction_flag = gs_input.split(',')[1] #write 'yes' to get the 2nd simpl
														#print('>>>>>____________________Computing bp simplification for',subfolder, funct, t, graph_type, min_,weighting_method_graph, weighting_method_simplification)
														i=0
														beta_discr = float(beta_discr)
														i+=1
														#updating_beta_discrete(beta_discr)

														Simp=graph_filtering_from_dat_files(subfolder,
																							t,
																							graph_type,
																							beta_discr,
																							funct,
																							min_,
																							btns_factor_source,
																							btns_factor_sink,
																							weighting_method_graph,
																							weighting_method_simplification,
																							source_sink[0],source_sink[1],
																							BP_weights,
																							reduction_flag)

														#print('simulation failed')
														#errors.append([funct,graph_type,weighting_method_graph,threshold,minimum,weighting_method_simplification,btns_factor])
														
														#os.chdir("../../Tests/")
						print('errors!',errors)
					


					
					############## GRAPH EXTRACTION (only for a single case) ######################
'''
					errors=[]
					for funct in ['flux']:
						for graph_type in ['3']: #','3']:
							for weighting_method_graph in ['ER']:
								if funct =='tdens':
									if 0<float(beta)<=1:
										t_list=[.5,.65]
									else:
										t_list=[0.001]
								else:
									if 0<float(beta)<=1:
										t_list=[.5,.65]
									else:
										t_list=[0.001]

								for threshold in t_list:
									subfolder='./runs/'+folder_name
									t=float(threshold)

									new_dir = subfolder+'/'+funct

									try:
										os.mkdir(new_dir)
									except OSError:
										print ("Creation of the directory %s failed." % new_dir)
									print('Flag and tolerance:',subfolder,t,graph_type,funct)
									G=getting_graphs(subfolder, t, graph_type,funct, weighting_method_graph) 

									for minimum in [0.001]:
										for weighting_method_simplification in ['IBP']:#, 'BPW', 'ER', None]:
											for btns_factor in [0.01]:#,0.05,0.1]:
												min_=float(minimum)
												btns_factor=float(btns_factor)
												#print('>>>>>____________________Computing bp simplification for',subfolder, funct, t, graph_type, min_,weighting_method_graph, weighting_method_simplification)
												try:
													Simp=BP_simplification(subfolder,t,graph_type,funct,min_,btns_factor,weighting_method_graph,weighting_method_simplification)
												except OSError: 
													print('simulation failed')
													errors.append([funct,graph_type,weighting_method_graph,threshold,minimum,weighting_method_simplification,btns_factor])
						print('errors!',errors)
					
					#print('This graph has:',str(len(G.nodes())),str(len(G[subfolder][str(t)][graph_type].edges())))
					#t_string= '%.0E' % decimal.Decimal( str(t) )
					#path_=subfolder+'/'+funct+'/'+funct+'_graph_t'+t_string+'_graph'+str(graph_type)+'.dat'
					#print('saving graph at',path_)
					#with open(path_, 'wb') as file:
					#    pkl.dump(G, file)

				##################FILE COPY####################

				#########BP SIMPLIFICATION#####

				#	for subfolder in ["./runs/" + folder_name]:
					
					funct='tdens'
					min_=float(minimum)
					graph_type= 1
					btns_factor=float(btns_factor)
					weighting_method = 'AVG'
					weighting_method_simplification = 'IBP'
					#subfolder = dst + "/" + folder_name
					print('computing bp simplification')
					Simp=BP_simplification(subfolder,t,graph_type,funct,min_,btns_factor,weighting_method,weighting_method_simplification)
'''						
##############

if __name__ == '__main__':
	new()

# todo: add parameter.txt prints