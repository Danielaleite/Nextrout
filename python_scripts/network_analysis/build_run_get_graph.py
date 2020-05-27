#!/usr/bin/env python


import click
import os
import shutil
import distutils
import numpy as np
import shutil
import distutils
from threading import Thread

@click.command()
@click.option('--refinements', prompt='What is the number of refinements?', help='Define the number of refinements.')
@click.option('--flag', prompt='What is your flag?', help='This is where you should put the name of the flag.')
#@click.option('--beta', prompt='Which value for beta?', help='This is where you should put your flag.')

def new(refinements,flag):#,inputs_crtl,folder_name):
	#def call_flag(refinements, flag):
	click.echo('Your inputs.ctrl file is ready! ')
	#click.echo('Your flag is %s!' %flag)
	file=open('assembly.ctrl','w+')

	file.write("# DOMAIN ####################################################### " + "\n")
	file.write('2d ! flag_domain (2d, 3d, surface)' + "\n")

	file.write('# MESH #########################################################' + "\n")
	file.write('rect_cnst path/grid.dat ! flag_grid' + "\n")
	file.write('16 ! ndiv' + "\n")
	file.write('%s' %refinements + ' ! nref' + "\n")

	file.write('# PREBUILD EXAMPLES ############################################' + "\n")
	file.write('1.5 extra_path ! flag_pflux ' + "\n")
	file.write('1.0 extra_path ! flag_pmass (gamma)' + "\n")
	file.write('1.0 extra_path ! flag_decay' + "\n")
	file.write('1.0 extra ! decay0' + "\n")
	file.write('rect_cnst  path/frog_source.dat ! flag_source' + "\n")
	file.write('rect_cnst  path/frog_sink2.dat ! flag_sink' + "\n")
	file.write('0 ! flag_normalize' + "\n")
	#file.write('%s' %flag + '!flag_tdens0 ' + "\n" )
	file.write('%s' %flag + ' ' + 'path/tdens0.dat !flag_tdens0 ' + "\n" )
	file.write('1.0 extra_path ! flag_kappa' + "\n")

	file.write('## TIME RANGE ###################################################' + "\n")
	file.write('0.0 ! tzero' + "\n")
	file.write('5.0e2 ! tmax (stop while loop)')
	#list = file.readlines()
	file.close()

	#shutil.copyfile('inputs.ctrl', 'Test_set/inputs.ctrl')

	folder_name = "%s" %flag
	def init(folder_name):
 		dir_otp="Test_set/"
 		os.chdir(dir_otp)
 		print(folder_name)
		new_dir="./runs/"+folder_name
		try:
			os.mkdir(new_dir)
		except OSError: 
			print ("Creation of the directory %s failed." % new_dir)

    
	command="./dmk_folder.py assembly " + './runs/'+folder_name + " assembly.ctrl "
	os.system(command)

	    #newdir2 = "./runs/"+folder_name.split("/")[-2]
	    #os.chdir('../')

	command="./dmk_folder.py run " + './runs/'+ folder_name  + " muffa.ctrl "
	os.system(command)

	command="./dmk_folder.py get-graph " +  './runs/'+ folder_name  + " 0.1 " + " 100000000000 "
	os.system(command)

	command="./dmk_folder.py vtk -a -tdens " +  './runs/'+ folder_name 
	os.system(command)

if __name__ == '__main__':
	new()