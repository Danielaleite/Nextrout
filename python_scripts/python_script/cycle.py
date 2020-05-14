#python

import numpy as np
import re
import os
from shutil import copyfile
import fileinput
import sys
sys.path.append('python_script/')
import python_script.common as common
import python_script.makejob as makejob
import itertools
import python_script.replace as rpl
import python_script.folder_structure as folder_module
from python_script.optvtk import optvtk


class change:
    def __init__(self,flag,file,column,values,label=None,labels=None):
        self.flag=str(flag)
        self.file=str(file)
        self.column=int(column)
        self.values=[str(value) for value in values]
        if (label==None) and (labels==None):
            self.labels=[
                str(self.flag)+str(value).replace('.','') 
                for value in self.values]
        else:
            if (label != None) :
                self.labels=[ 
                    str(label) +
                    str(value).replace('.','') 
                    for value in self.values]
            if (labels != None):
                self.labels=labels

#
# get locations
#
folder_runs=os.path.abspath('./runs/')

file_ctrl='inputs.ctrl'

file_location = 'location_repository.txt'
flocation = open(file_location, "r")
lines = flocation.readlines()
repo_path=str(lines[0]).rstrip()
flocation.close() 


def cycle(changes,order,folder_base=None):
    combination=list(
        itertools.product(*[range(len(change.values)) for change in changes]))

    for comb in combination:
        #
        # define folder name
        #      
        values=[change.values[comb[j]] for j, change in enumerate(changes)]
        labels=[change.labels[comb[j]] for j, change in enumerate(changes)]
        # print values
        # print labels
        # print "_".join(labels)
        
        foldername=folder_base + "_".join(labels)
        folderpath=str(folder_runs)+'/'+str(foldername)

        #
        # define subfolder structrure
        # 
        folder = folder_module.folder_structure(folderpath)

        if ( order=='test'):
            print ' '
            print 'LABELS     = ',labels
            print 'VALUES     = ',values
            print 'FOLDER_NAME= ',folder.name

        if ( order=='inputs'):
            print ' '
            print 'LABELS     = ',labels
            print 'VALUES     = ',values
            print 'FOLDER_NAME= ',folder.name
        
            #
            # replace values
            #
            for ich in range(len(changes)):
                rpl.replace(changes[ich].flag, 
                            changes[ich].file, 
                            changes[ich].column, 
                            changes[ich].values[comb[ich]])

                
            #
            # assembly inputs
            # 
            assembly_inputs(foldername)

        if ( order=='folders'):
            #
            # create folder and subfolders
            # 
            folder.mkdirs()        
        
        elif ( order=='run'):
            #
            # copy controls  file 
            #
            command='cp muffa.fnames '+ str(folder.name)
            os.system(command)
            command='cp muffa.ctrl '  + str(folder.input)
            os.system(command)

            #
            # mov into the folder and launch
            #
            os.chdir(str(folder.name))
            command='../../../repositories/muffe_p1p0/code/muffa.out'
            os.system(command)
            os.chdir('../../launch_controls')

        
        elif ( order=='timevtk'):
            timevtk(str(folder.name))
        
        elif ( order=='optvtk'):
            optvtk([folder.name])

        elif ( order=='makejob'):
            makejob.makejob(str(os.path.abspath('jobs/base.job')),
                            str(os.getcwd()),
                            foldername, folder.name,&
                            'short')

        elif ( order=='figures'):
            #
            # define figures path and mkdir if necessary
            # 
            folder_figures_path  = '../runs/figures/'
            if not os.path.exists( folder_figures_path ):
                os.makedirs( folder_figures_path )
                
            input = open('base.session')
            output = open('work.session','w')

            for line in input:
                output.write(re.sub(r'folder_name[^\/]*/', str(folder_path), line))

            input.close()
            output.close()
        
            #os.mkdir(str(folder_figures_path))
            command=('visit -cli -nowin -s ./pyhton_script/last.py')
            os.system(command)


            print 'ninputs=', len(combination)


#
# list of changes defined by
# 1 - flag mathcing on file
# 2 - file name
# 3 - index of column (start-from-0 enumeration)
# 4 - list of values reqires
#

order=sys.argv[1]
changes=[]
changes.append(change('flag_domain','inputs.ctrl',0,['2d','3d'],'setting'))
cycle(changes,order,'prova_')

#changes.append(change('pflux0','inputs.ctrl',0,['1.50','2.00','3.00'],'pflux'))
# changes.append(change('flag_tdens0','inputs.ctrl',0,['one','par0505','ybranch'],'tdini'))
# changes.append(change('ndiv','inputs.ctrl',0,[63],'ndiv'))

# article setting
#changes.append(change('pflux0','inputs.ctrl',0,['1.30','1.70','1.80','1.90'],'pflux'))
#changes.append(change('flag_tdens0','inputs.ctrl',0,['one','par0505'],'tdini'))
#changes.append(change('ndiv','inputs.ctrl',0,[63,64],'ndiv'))
#


#changes.append(change('alpha','info_grids/y_branch.dat',0,[0.0,0.9],'alpha'))
#changes.append(change('pflux0','inputs.ctrl',0,['1.30','1.70','1.80','1.90'],'pflux'))

#changes.append(change('pflux','inputs.ctrl',0,[
#'1.05','1.10','1.20','1.30','1.40','1.50',
#   '1.60','1.70'],#,'1.20','1.30','1.40','1.50'],
#  label='pflux'))

# changes.append(change('flag_tdens0','inputs.ctrl',1,[
#     '../bergen/tdens0_2d.dat',
#     '../bergen/tdens0_discon_2dblur.dat',
#     '../bergen/tdens0_discon_2d.dat'],
#     labels=[
#         'td0',
#         'td0discon',
#         'td0disconblur']))




#     '1.40','1.50','1.60','1.70',
#     '1.80','1.90','2.00','2.30',

#changes.append(change('ndiv','inputs.ctrl',0,[63,64],'ndiv'))




#
# define optional root of folder name
# folder=['roots']+listof('flag'+'value')
#

 
#labels=[change.label for change in changes]
#print labels
#print len(changes[-1].label)

        
        
        
    
    
