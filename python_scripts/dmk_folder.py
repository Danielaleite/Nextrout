#!/usr/bin/env python
import click
import os
import fnmatch
import glob

import python_script.dmk_folder_structure as dmk 
import python_script.utilities as util
from shutil import copyfile

#
# include path to muffe_p1p0
#
file = open('location_repository.txt', 'r')
repo = file.read().strip()+'/'
file.close()
import sys
sys.path.append(str(repo))
#sys.path.append(str(repo)+'/preprocess')
sys.path.append(str(repo)+'../globals/python_timedata')
import timedata as timedata 


def find_path_file(abspath,flag_input):
    fnames=abspath+'/muffa.fnames'
    path=util.search_read(fnames,flag_input,0)
    if ( path[0]=='/') :
        return path
    else :
        return os.path.normpath(abspath+'/'+path)


@click.group()
def cli():
    pass

@cli.command()
@click.argument('fname')
def init(fname):
    """Create empty folder"""
    folder=dmk.dmk_folder_structure(fname)
    folder.mkdirs()
    click.echo('Initialized dmk folder %s' % fname)


@cli.command()
@click.argument('path_original')
@click.argument('path_copy')
def cp(path_original, path_copy):
    """ Create a copy with muffa.fnames pointing 
    \b
    to the inputs of the orignal
    """
    #
    # init. original and copy structure (creates copy folders if required )
    #
    original=dmk.dmk_folder_structure(path_original)
    copy=dmk.dmk_folder_structure(path_copy)
    if ( not os.path.exists(copy.abspath) ):
        click.echo('making folders')
        copy.mkdirs()
    
    #
    # get relative path
    #
    rel_path=(os.path.relpath(original.abspath, copy.abspath))
    
    #
    # copy fnames files
    #
    src=original.abspath+'/muffa.fnames'
    dst=copy.abspath+'/muffa.fnames'
    copyfile(src,dst)
    
    #
    # reassign all input file via muffa.fnames
    #
    input_flags=['flag_grid',
                 'flag_subgrid',
                 'flag_parent',
                 'flag_kappa',
                 'flag_pflux',
                 'flag_pmass',
                 'flag_decay',
                 'flag_tdens0',
                 'flag_rhs_integrated',
                 'flag_optdens',
                 'flag_penalty_function',
                 'flag_penalty_weight',
                 'flag_penalty_factor',
                 'flag_dirichlet']
    
    #
    # read directory
    #
    directory=util.search_read(str(dst), 'flag_directory', 0)

    #
    # read file path
    #
    for flag in input_flags:
        print(flag)
        original_path=util.search_read(str(dst), flag, 0)
        if (original_path != '') :
            click.echo(flag +' '+ original.abspath+'/'+original_path )
            if ( os.path.exists(original.abspath+'/'+original_path) ):
                if (original_path[0] == '/'):
                    #
                    # absolute path case
                    #
                    redirect=original_path
                else:
                    #
                    # relarive path
                    #
                    redirect=rel_path+'/'+original_path
                    util.replace(str(dst),flag,0,redirect)

    click.echo('Folder %s copy into %s' % (original.path,copy.path) )


@cli.command()
@click.argument('folder_name')
@click.option('flag_option', '-o', flag_value='o',
              default=True,
              help='First  Option: Build  Optimal VTK files')
@click.option('flag_option', '-a', flag_value='a',
              help='First  Option: All time VTK files')

@click.option('flag_data', '-all', flag_value='all',
              default=True,
              help='Second Option: Build Tdens and Pot VTK files')
@click.option('flag_data', '-tdens', flag_value='tdens',
              help='Second Option: Build only Tdens VTK files')
@click.option('flag_data', '-pot', flag_value='pot',
              help='Second Option: Build only Pot VTK files')
@click.option('flag_data', '-gradavg', flag_value='gradavg',
              help='Second Option: Build only Pot VTK files')


def vtk(folder_name,flag_option,flag_data):
    """Make vtk files of optimal data"""
    #
    # check existence
    #
    if ( not os.path.exists(folder_name) ) :
        clik.echo( 'folder < ' + str(folder_name) + ' > not found')
        return
    #
    # program used
    #
    file = open('location_repository.txt', 'r')
    repo = file.read().strip()+'/'
    file.close()
    program=repo+'../geometry/timedata2vtk/timedata2vtk.out '
    program=os.path.normpath(program)
    #
    # get input files
    #
    folder=dmk.dmk_folder_structure(folder_name)
    grid_path=find_path_file(folder.abspath,'flag_grid')
    subgrid_path=find_path_file(folder.abspath,'flag_subgrid')
   
    if ( flag_option == 'o'):
        if ( flag_data == 'tdens') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(grid_path)    + ' ' +
                     str(folder.result)+'/opt_tdens.dat' + ' ' +
                     str(folder.output_vtk)+'/')
            os.system(command)
        
        if ( flag_data == 'pot') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(subgrid_path)    + ' ' +
                     str(folder.result)+'/opt_pot.dat' + ' ' +
                     str(folder.output_vtk)+'/')
            os.system(command)

        if ( flag_data == 'gradavg') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(grid_path)    + ' ' +
                     str(folder.result)+'/opt_nrm_grad_avg.dat' + ' ' +
                     str(folder.output_vtk)+'/')
            os.system(command)
        if ( flag_data == 'velprj') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(grid_path)    + ' ' +
                     str(folder.result)+'/opt_vel_prj.dat' + ' ' +
                     str(folder.output_vtk)+'/')
            os.system(command)

    if ( flag_option == 'a'):
        if ( flag_data == 'tdens') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(grid_path)    + ' ' +
                     str(folder.result)+'/tdens.dat' + ' ' +
                     str(folder.output_vtk+'/'))
            os.system(command)
        
        if ( flag_data == 'pot') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(subgrid_path)    + ' ' +
                     str(folder.result)+'/pot.dat' + ' ' +
                     str(folder.output_vtk+'/'))
            os.system(command)
        if ( flag_data == 'gradavg') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(grid_path)    + ' ' +
                     str(folder.result)+'/nrm_grad_avg.dat' + ' ' +
                     str(folder.output_vtk)+'/')
            os.system(command)
        if ( flag_data == 'velprj') or ( flag_data == 'all'): 
            command=(str(program)      + ' ' +
                     str(grid_path)    + ' ' +
                     str(folder.result)+'/vel_prj.dat' + ' ' +
                     str(folder.output_vtk)+'/')
            os.system(command)
        

@cli.command()
@click.argument('folder_name')
@click.option('flag_option', '-o', flag_value='o',
              default=True,
              help='Option 1 : Build Optimal Extra Data (Default)')
@click.option('flag_option', '-a', flag_value='a',
              help='Option 2 : Build Evolution Extra Data')

# @click.option('flag_data', '-all', flag_value='all',
#               default=True,
#               help='Second Option: Build all extradata')
# @click.option('flag_data', '-grad', flag_value='grad',
#               help='Second Option: Build gradient of pot')
# @click.option('flag_data', '-grad_prj', flag_value='grad_prj',
#               help='Second Option: Build gradient of pot projected of tdens grid')
# @click.option('flag_data', '-vel', flag_value='vel',
#               help='Second Option: Build vel')
def extradata(folder_name,flag_option):
    """Build gradient, projected gradient, vel"""     
    folder=dmk.dmk_folder_structure(folder_name)


    #
    # define paths to input file (search in muffa.fnames) 
    #
    grid_path=find_path_file(folder.abspath,'flag_grid')
    subgrid_path=find_path_file(folder.abspath,'flag_subgrid')
    parent_path=find_path_file(folder.abspath,'flag_parent')
   
    if ( flag_option == 'o'):
        #
        # compute graddient of p1
        #
        program=repo+'../p1galerkin/eval_gradient/eval_gradient.out'
        command=(str(program)      + ' ' +
                 str(subgrid_path) + ' ' +
                 str(folder.result)+'/opt_pot.dat' + ' ' +
                 str(folder.result)+'/opt_grad.dat')
        os.system(command)
        click.echo('Optimal Gradient           - Done')

        #
        # compute projected gradient
        #
        
        # average gradient of grid tdens
        program=repo+'../geometry/project_timedata/project_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/opt_grad.dat'+ ' ' +
                 str(subgrid_path) + ' ' +
                 str(grid_path)    + ' ' +
                 str(parent_path)  + ' ' +
                 str(folder.result)+'/opt_grad_prj.dat')
        os.system(command)
        click.echo('Optimal Gradient Projected - Done')

        #click.echo('Optimal Projected Tdens        - Done')
        # multiply tdens_prj * gradient
        # interpolate opt_tdens
        program=os.path.abspath(repo+
                                '../geometry/interpolate_timedata/interpolate_timedata.out ')
        command=(str(program)      + ' ' +
                 subgrid_path      + ' ' + 
                 parent_path      + ' ' +
                 str(folder.result)+'/opt_tdens.dat'+' ' +
                 str(folder.result)+'/opt_tdens_inter.dat')
        os.system(command)
        program=repo+'../globals/multiply_timedata/multiply_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/opt_tdens_inter.dat'+' ' +
                 str(folder.result)+'/opt_grad.dat' +' ' +
                 str(folder.result)+'/temp.dat')
        os.system(command)
        # Reverse Vel
        program=repo+'../globals/scale_timedata/scale_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/temp.dat'+' ' +
                 ' -1.0  ' +
                 str(folder.result)+'/opt_vel.dat')
        os.system(command)
        click.echo('Optimal Vel                - Done')
        #os.remove(str(folder.result)+'/opt_tdens_inter.dat')
                
        
        # multiply tdens * prjected gradient
        program=repo+'../globals/multiply_timedata/multiply_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/opt_tdens.dat'+' ' +
                 str(folder.result)+'/opt_grad_prj.dat' +' ' +
                 str(folder.result)+'/temp.dat')
        os.system(command)
        program=repo+'../globals/scale_timedata/scale_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/temp.dat'+' ' +
                 ' -1.0  ' +
                 str(folder.result)+'/opt_vel_prj.dat')
        os.system(command)
        click.echo('Optimal Vel Projected      - Done')
    
        #
        # remove work file
        #
        #os.remove(str(folder.result)+'/temp.dat')
        

    if ( flag_option == 'a'):
        #
        # compute graddient of p1
        #
        program=repo+'../p1galerkin/eval_gradient/eval_gradient.out'
        command=(str(program)      + ' ' +
                 str(subgrid_path) + ' ' +
                 str(folder.result)+'/pot.dat' + ' ' +
                 str(folder.result)+'/grad.dat')
        os.system(command)
        click.echo('Gradient           - Done')

        #
        # compute projected gradient
        #
        
        # average gradient of grid tdens
        program=repo+'../geometry/project_timedata/project_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/grad.dat'+ ' ' +
                 str(subgrid_path) + ' ' +
                 str(grid_path)    + ' ' +
                 str(parent_path)  + ' ' +
                 str(folder.result)+'/grad_prj.dat')
        os.system(command)
        click.echo('Gradient Projected - Done')

        #click.echo('Optimal Projected Tdens        - Done')
        # multiply tdens_prj * gradient
        # interpolate opt_tdens
        program=os.path.abspath(repo+
                                '../geometry/interpolate_timedata/interpolate_timedata.out ')
        command=(str(program)      + ' ' +
                 subgrid_path      + ' ' + 
                 parent_path      + ' ' +
                 str(folder.result)+'/tdens.dat'+' ' +
                 str(folder.result)+'/tdens_inter.dat')
        os.system(command)
        program=repo+'../globals/multiply_timedata/multiply_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/tdens_inter.dat'+' ' +
                 str(folder.result)+'/grad.dat' +' ' +
                 str(folder.result)+'/temp.dat')
        os.system(command)
        # Reverse Vel
        program=repo+'../globals/scale_timedata/scale_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/temp.dat'+' ' +
                 ' -1.0  ' +
                 str(folder.result)+'/vel.dat')
        os.system(command)
        click.echo('Vel                - Done')
                
        
        # multiply tdens * prjected gradient
        program=repo+'../globals/multiply_timedata/multiply_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/tdens.dat'+' ' +
                 str(folder.result)+'/grad_prj.dat' +' ' +
                 str(folder.result)+'/temp.dat')
        os.system(command)
        program=repo+'../globals/scale_timedata/scale_timedata.out '
        command=(str(program)      + ' ' +
                 str(folder.result)+'/temp.dat'+' ' +
                 ' -1.0  ' +
                 str(folder.result)+'/vel_prj.dat')
        os.system(command)
        click.echo('Vel Projected      - Done')
    
        #
        # remove work file
        #
        os.remove(str(folder.result)+'/temp.dat')


@cli.command()
@click.argument('folder_name')
@click.argument('flag_input')
@click.argument('input_file_path')
def assign(folder_name, flag_input, input_file_path):
    """Assign an input file"""
    #
    # get relative path of file with respect to folder
    #
    folder=dmk.dmk_folder_structure(folder_name)
    rel_path=os.path.relpath(str(input_file_path), str(folder.abspath))
    util.replace(folder.abspath+'/muffa.fnames',flag_input,0,rel_path)


@cli.command()
@click.argument('folder_name')
def optfield(folder_name):
    """Compute optimal field (Evans-Gangbo-1999)"""
    #
    # path to program used
    #
    program=os.path.normpath(repo+
                             '../build_optimal_vector_field/code/optimal_vector_field.out')
    
    #
    # get input files
    #
    folder=dmk.dmk_folder_structure(folder_name)
    grid_path=find_file_path(folder.abspath,'flag_grid')
    subgrid_path=find_file_path(folder.abspath,'flag_subgrid')
    parent_path=find_file_path(folder.abspath,'flag_parent')
    command=(program     + ' ' +
             grid_path   + ' ' + 
             subgrid_path+ ' ' +
             parent_path+ ' ' +
             str(paths.input) +'forcing.dat'+ ' ' +
             str(paths.result)+'opt_tdens.dat'+ ' ' +
             str(paths.result)+'opt_pot.dat' + ' ' +
             str(paths.result)+'opt_field.dat')
    os.system(command)
    
@cli.command()
@click.argument('folder_name')
@click.argument('file_ctrl')
def assembly(folder_name,file_ctrl):
    """Use muffe_p1p0/preprocess examples"""
    folder=dmk.dmk_folder_structure(folder_name)
    folder.mkdirs()
    # repository location
    file = open('location_repository.txt', 'r')
    repo = os.path.abspath(file.read().strip())
    file.close()

    #
    # copy folder_names file 
    #
    src=repo+'/code/muffa.fnames'
    dst=folder.abspath+'/muffa.fnames'
    copyfile(src, dst)
    
    #
    # get assembler location
    # 
    absctrl=os.path.abspath(file_ctrl)
    domain=util.search_read(absctrl,'flag_domain',0)
    assembly_folder = repo+'/preprocess/'+domain+'_assembly/'
    
    #
    # mov into assembly folder and create inputs
    #
    cwd=os.getcwd()
    os.chdir(assembly_folder)
    print(assembly_folder)
    command=('python main.py'  +' '+
             str(absctrl) +' '+
             str(folder.input)  +' '+
             str(folder.input_vtk) )
    click.echo( command)
    os.system(command)
    os.chdir(cwd)

@cli.command()
@click.argument('folder_name')
@click.argument('ctrlfile')
def run(folder_name,ctrlfile):
    """Run simulation given folder and control file"""
    folder=dmk.dmk_folder_structure(folder_name)
    
    # repository location
    file = open('location_repository.txt', 'r')
    repo = os.path.abspath(file.read().strip())
    file.close()
    program=repo+'/code/muffa.out'

    #
    # clean possibly existing data
    # 
    for fl in glob.glob(folder.result+'/*.dat'):
        os.remove(fl)
    for fl in glob.glob(folder.timefun+'/*.dat'):
        os.remove(fl)
    for fl in glob.glob(folder.output_vtk+'/*.vtk'):
        os.remove(fl)
        
        


    #
    # copy control file
    #
    dst=folder.input+'/muffa.ctrl'
    copyfile(os.path.abspath(ctrlfile),dst )



    #
    # mov into assembly folder and create inputs
    #
    cwd=os.getcwd()
    os.chdir(folder.abspath)
    print(os.getcwd())
    os.system(program)
    os.chdir(cwd)

@cli.command()
@click.argument('folder_name')
def err_distance(folder_name):
    """Compute err=pot-exact_pot (after centering)"""
    #
    # get relative path of file with respect to folder
    #
    folder=dmk.dmk_folder_structure(folder_name)
    
    # read optimal solution and get point dist(p)=0
    opt_pot=timedata.read_steady_timedata(str(folder.input+'/opt_pot.dat'))
    clik.echo( opt_pot.shape)
    for i in range(len(opt_pot)):
        if (opt_pot[i] == 0.0 ):
            inode_sub = i
            break

    clik.echo('inode='+str(inode_sub))

    # find the corresponding value
    opt_pot=timedata.read_steady_timedata(str(folder.result+'/opt_pot.dat'))
    shift=opt_pot[inode_sub]
    
    clik.echo( shift)

    #
    # shift optimal potential (otp_pot_shifted.dat) 
    # and optimal potential (pot_shifted.dat)
    # and compute errors
    
    #
    # program used
    #
    program=os.path.normpath(repo+
                             '../globals/axpy_timedata/axpy.out')
    #
    # shift potentials
    #
    command=(program + ' ' +
             '-1.0 '+str(shift)+ ' ' +
             str(folder.result)+'/opt_pot.dat' + ' ' +
             str(folder.result)+'/opt_pot_shifted.dat')
    clik.echo(command)
    os.system(command)

    # command=(program + ' ' +
    #          '-1.0 '+str(shift)+ ' ' +
    #          str(folder.result)+'/pot.dat' + ' ' +
    #          str(folder.result)+'/pot_shifted.dat')
    # os.system(command)
      
    command=(program + ' ' +
             '-1.0 ' +
             str(folder.input+'/opt_pot.dat') + ' ' +
             str(folder.result+'/opt_pot_shifted.dat')+ ' ' +
             str(folder.result+'/err_opt_pot_shifted.dat'))
    os.system(command)

    # command=(program + ' ' +
    #          '-1.0 '+
    #          str(folder.input+'/opt_pot.dat') + ' ' +
    #          str(folder.result+'/pot_shifted.dat')+ ' ' +
    #          str(folder.result)+'/err_pot_shifted.dat')
    # os.system(command)

    #
    # program used
    #
    program=os.path.normpath(repo+
                             '../geometry/lp_norm_timedata/lp_norm.out')
    
    command=(program + ' ' +
             str(folder.input)+'/subgrid.dat' + ' ' +
             str(folder.result)+'/err_opt_pot_shifted.dat' + ' ' +
             ' 2.0')
    os.system(command)

# @cli.command()
# @click.argument('folder_name')
# @click.option('flag_option', '-tdens', flag_value='tdens',
#               default=True,
#               help='Compute errors w.r.t to Optimal Transport Density')
# @click.option('flag_option', '-pot', flag_value='pot',
#               help='First Compute errors w.r.t to Optimal Transport Density')   


# def error(folder_name):
#     """Compute err=pot-exact_pot (after centering)"""
#     #
#     # get relative path of file with respect to folder
#     #
#     folder=dmk.dmk_folder_structure(folder_name)
    
#     # read optimal solution and get point dist(p)=0
#     opt_pot=timedata.read_steady_timedata(str(folder.input+'/opt_pot.dat'))
#     clik.echo( opt_pot.shape)
#     for i in range(len(opt_pot)):
#         if (opt_pot[i] == 0.0 ):
#             inode_sub = i
#             break

#     clik.echo('inode='+str(inode_sub))

#     # find the corresponding value
#     opt_pot=timedata.read_steady_timedata(str(folder.result+'/opt_pot.dat'))
#     shift=opt_pot[inode_sub]
    
#     clik.echo( shift)

#     #
#     # shift optimal potential (otp_pot_shifted.dat) 
#     # and optimal potential (pot_shifted.dat)
#     # and compute errors
    
#     #
#     # program used
#     #
#     program=os.path.normpath(repo+
#                              '../globals/axpy_timedata/axpy.out')
#     #
#     # shift potentials
#     #
#     command=(program + ' ' +
#              '-1.0 '+str(shift)+ ' ' +
#              str(folder.result)+'/opt_pot.dat' + ' ' +
#              str(folder.result)+'/opt_pot_shifted.dat')
#     clik.echo(command)
#     os.system(command)

#     # command=(program + ' ' +
#     #          '-1.0 '+str(shift)+ ' ' +
#     #          str(folder.result)+'/pot.dat' + ' ' +
#     #          str(folder.result)+'/pot_shifted.dat')
#     # os.system(command)
      
#     command=(program + ' ' +
#              '-1.0 ' +
#              str(folder.input+'/opt_pot.dat') + ' ' +
#              str(folder.result+'/opt_pot_shifted.dat')+ ' ' +
#              str(folder.result+'/err_opt_pot_shifted.dat'))
#     os.system(command)

#     # command=(program + ' ' +
#     #          '-1.0 '+
#     #          str(folder.input+'/opt_pot.dat') + ' ' +
#     #          str(folder.result+'/pot_shifted.dat')+ ' ' +
#     #          str(folder.result)+'/err_pot_shifted.dat')
#     # os.system(command)

#     #
#     # program used
#     #
#     program=os.path.normpath(repo+
#                              '../geometry/lp_norm_timedata/lp_norm.out')
    
#     command=(program + ' ' +
#              str(folder.input)+'/subgrid.dat' + ' ' +
#              str(folder.result)+'/err_opt_pot_shifted.dat' + ' ' +
#              ' 2.0')
#     os.system(command)



@cli.command()
@click.argument('folder_name')
@click.argument('min')
@click.argument('max')
@click.option('n', '--o',flag_value='o',type=str,
              default=True,
              help='Build optimal of graph')
@click.option('n', '--a', flag_value='a',help='Build sequence of graphs')


def get_graph(folder_name,min,max,n):
    """Build graph thresholding opt_tdens or thresholding all the obtained tdens"""
    folder=dmk.dmk_folder_structure(folder_name)
    if ( n =='a'):
        # This divides the tdens file in separated files, one for each iter number
        ### https://stackoverflow.com/questions/546508/how-can-i-split-a-file-in-python
        input = open(folder.result+'/opt_tdens.dat', 'r').read().split('\n')
        firstLine=input[0]
        dim_=int(input[0].split()[1])
        outputBase = folder.result+'/tdens'

        # This is shorthand and not friendly with memory
        # on very large files (Sean Cavanagh), but it works.
        input = open(folder.result+'/tdens.dat', 'r').read().split('\n')

        splitLen = dim_+2
        at = 0
        for lines in range(0, len(input), splitLen):
            # First, get the list slice
            outputData = input[lines+1:lines+splitLen+1]
            # Now open the output file, join the new slice with newlines
            # and write it out. Then close the file.
            output = open(outputBase + str(at) + '.dat', 'w')
            output.write(firstLine+'\n')
            output.write('\n'.join(outputData))
            output.write('\n time     1e+30     ')
            output.close()
            # Increment the counter
            at += 1
    
    #
    # build graph of cell connection
    #
    full_graph=str(folder.input+'/graph_cell.dat')

    if ( not os.path.exists(full_graph) ):
        program=os.path.normpath(repo+
                                 '../geometry/grid2graph/grid2graph.out')
        command=(program + ' ' +
                 str(folder.input+'/grid.dat')+ ' ' +
                 full_graph + ' cell')
        os.system(command)

    if n == 'a':
        # -a stands for 'all'
        clik.echo( '-a. Building all the graphs')
        program1=os.path.normpath(repo+
                                 '../geometry/selection_grid/selection_grid.out')
        program2=os.path.normpath(repo+'../geometry/grid2vtk/grid2vtk.out ')
        for i in range(at-1):
            
            command=(program1 + ' ' +
                 full_graph + ' ' +
                 str(folder.result+'/tdens'+str(i)+'.dat') + ' ' +
                 str(min) + ' ' + str(max)+ ' ' +
                 str(folder.result+'/graph'+str(i)+'.dat') + ' ' +
                 str(folder.result+'/selector_graph.dat'))
    
            os.system(command)
            #
            # print vtk
            #
            
            command=(str(program2)      + ' ' +
                     str(folder.result+'/graph'+str(i)+'.dat') + ' ' +
                     str(folder.output_vtk)+'/graph'+str(i)+'.vtk')
            os.system(command)

    elif n=='o':
        # -o stand for 'optimal'
        click.echo( '-o. Building optimal graphs')
    
        #
    # build graph thresholding opt_tdens
    #
        program=os.path.normpath(repo+
                                 '../geometry/selection_grid/selection_grid.out')
        command=(program + ' ' +
                 full_graph + ' ' +
                 str(folder.result+'/opt_tdens.dat') + ' ' +
                 str(min) + ' ' + str(max)+ ' ' +
                 str(folder.result+'/opt_graph.dat') + ' ' +
                 str(folder.result+'/selector_graph.dat'))
    
        os.system(command)

        #
        # print vtk
        #
        program=os.path.normpath(repo+'../geometry/grid2vtk/grid2vtk.out ')
        command=(str(program)      + ' ' +
                 str(folder.result+'/opt_graph.dat') + ' ' +
                 str(folder.output_vtk)+'/opt_graph.vtk')
        os.system(command)
    else:
        clik.echo('-'+str(n)+'. Not a valid input.')


@cli.command()
@click.argument('folder_name')
@click.argument('selector_data')
@click.argument('min')
@click.argument('max')
@click.argument('folder_selection')

def select(folder_name,selector_data,min,max,folder_selection):
    """ Define new folder selecting portion of data"""

    #
    # create empty folder
    # 
    folder=dmk.dmk_folder_structure(folder_name)
    folder_out=dmk.dmk_folder_structure(folder_selection)
    folder_out.mkdirs()
    src=repo+'code/muffa.fnames'
    dst=folder_out.abspath+'/muffa.fnames'
    click.echo(src)
    click.echo(dst)
    copyfile(src, dst) 


    #
    # get path of spatial discretization
    #
    grid_path=find_path_file(folder.abspath,'flag_grid')    
    subgrid_path=find_path_file(folder.abspath,'flag_subgrid')    
    parent_path=find_path_file(folder.abspath,'flag_parent')    
    click.echo(parent_path)
    
    #
    # copy selector file
    #
    dst=folder_out.input+'/selector_data_grid.dat'
    copyfile(selector_data, dst) 

    #
    # project selector to subgrid
    #
    program=os.path.normpath(repo+
                             '../geometry/interpolate_timedata/interpolate_timedata.out')
    command=(program   + ' ' +
             subgrid_path + ' ' +
             parent_path + ' ' +
             str(selector_data) + ' ' +
             str(folder_out.input+'/selector_data_subgrid.dat'))
    os.system(command)
    
    #
    # build project grid and get selctor file
    #
    program=os.path.normpath(repo+
                             '../geometry/selection_grid/selection_grid.out')
    command=(program   + ' ' +
             grid_path + ' ' +
             str(selector_data) + ' ' +
             str(min) + ' ' + str(max)+ ' ' +
             str(folder_out.input+'/grid.dat') + ' ' +
             str(folder_out.input+'/selector_cell_grid.dat') + ' ' +
             str(folder_out.input+'/selector_node_grid.dat'))
    os.system(command)
                            
    command=(program   + ' ' +
             subgrid_path + ' ' +
             str(folder_out.input+'/selector_data_subgrid.dat') + ' ' +
             str(min) + ' ' + str(max)+ ' ' +
             str(folder_out.input+'/subgrid.dat') + ' ' +
             str(folder_out.input+'/selector_cell_subgrid.dat') + ' ' +
             str(folder_out.input+'/selector_node_subgrid.dat'))
    os.system(command)

    program=os.path.normpath(repo+
                             '../geometry/selection_parent/selection_parent.out')
    
    #
    # buidl new parent file
    #
    parent_path=find_path_file(folder.abspath,'flag_parent')
    path_selected=find_path_file(folder_out.abspath,'flag_parent')
    click.echo( path_selected)
    command=(program   + ' ' +
             parent_path + ' ' +
             str(folder_out.input+'/selector_node_grid.dat') + ' ' +
             str(folder_out.input+'/selector_cell_grid.dat') + ' ' +
             str(folder_out.input+'/selector_node_subgrid.dat') + ' ' +
             str(folder_out.input+'/selector_cell_subgrid.dat') + ' ' +
             path_selected)
    os.system(command)


    
    #
    # read file with scalar quantities and just copy
    #
    program=os.path.normpath(repo+
                             '../globals/selection_timedata/selection_timedata.out')
    input_flags=[ 'flag_pflux',
                  'flag_pmass',
                  'flag_decay',
                  'flag_penalty_factor']
    for flag in input_flags:
        path_original=find_path_file(folder.abspath,flag)
        path_selected=find_path_file(folder_out.abspath,flag)
        
        #click.echo(path_original)
        #click.echo(path_selected)
        if ( os.path.exists(path_original) ):
            copyfile(path_original, path_selected)
    
    #
    # read file with grid quantities and select
    #              
    input_flags=['flag_kappa',
                 'flag_tdens0',
                 'flag_optdens',
                 'flag_penalty_function',
                 'flag_penalty_weight']
    for flag in input_flags:
        path_original=find_path_file(folder.abspath,flag)
        path_selected=find_path_file(folder_out.abspath,flag)
        #click.echo(path_original)
        #click.echo(path_selected)
        if ( os.path.exists(path_original) ):
            command=(program   + ' ' +
                     path_original + ' ' +
                     str(folder_out.input+'/selector_cell_grid.dat')+ ' ' +
                     path_selected)
            os.system(command)
    input_flags=['flag_rhs_grid_integrated']
    for flag in input_flags:
        path_original=find_path_file(folder.abspath,flag)
        path_selected=find_path_file(folder_out.abspath,flag)
        #click.echo(path_original)
        #click.echo(path_selected)
        if ( os.path.exists(path_original) ):
            command=(program   + ' ' +
                     path_original + ' ' +
                     str(folder_out.input+'/selector_node_grid.dat')+ " " +
                     path_selected)
            click.echo(command)
            os.system(command)

                 

                             
    input_flags=[ 'flag_rhs_integrated',
                  'flag_dirichlet']
    for flag in input_flags:
        path_original=find_path_file(folder.abspath,flag)
        path_selected=find_path_file(folder_out.abspath,flag)
        #click.echo(path_original)
        #click.echo(path_selected)
        if ( os.path.exists(path_original) ):
            command=(program   + ' ' +
                     path_original + ' ' +
                     str(folder_out.input+'/selector_node_subgrid.dat')+ " " +
                     path_selected)
            os.system(command)

@cli.command()
@click.argument('folder_name')
@click.option('flag_option', '-t', flag_value='t',
              default=True,
              help='First  Option: Time vs Data')
@click.option('flag_option', '-a', flag_value='a',
              help='First  Option: Alogorithm time vs Data')

def timevsdata(folder_name,flag_option):
    """Creates file with x_column=time y_column=data"""

    folder=dmk.dmk_folder_structure(folder_name)
    if ( flag_option == 't') :
        time_file=folder.timefun+'/time.dat'
    else:
        time_file=folder.timefun+'/algorithm.dat'
    
    #
    # select data file
    #
    timefun_files=os.listdir(os.path.abspath(folder.timefun))
    timefun_files=[ x for x in timefun_files if 'time' not in x]
    timefun_files=[ x for x in timefun_files if 'algorithm' not in x]


    for f in timefun_files:
        if ( flag_option == 't') :
            file_time_vs_data=(folder.timefun +
                               '/time_vs_'+os.path.basename(f))
        else:
            file_time_vs_data=(folder.timefun +
                               '/algorithm_vs_'+os.path.basename(f))
        data_file=os.path.abspath(folder.timefun)+'/'+f
        util.print_time(time_file, data_file, file_time_vs_data)

@cli.command()
@click.argument('folder_name')
@click.argument('value')
                #help='Real value required')
                
@click.argument('file_out')
                #help='Path for required data')
@click.option('flag_option', '-v', flag_value='v',
              default=True,
              help='First  Option: Get data at given variation of tdens ')
@click.option('flag_option', '-t', flag_value='t',
              help='First  Option: Get data at given time ')

def get_data(folder_name,flag_option,value,file_out):
    """Get tdens at given time"""
        
    folder=dmk.dmk_folder_structure(folder_name)

    value=float(value)

    if (flag_option == 't') :
        time=value
    if (flag_option == 'v') :
        #
        # find line in var_tdens
        #
        filepath=folder.timefun+'/var_tdens.dat'
        filein=open(str(filepath), 'r')
        input_lines = filein.readlines()
        ndata=len(input_lines)
        found=False
        i=0
        while ((not found ) and (i < ndata) ) :
            var=float(input_lines[i].rstrip())
            if ( var < value) :
                found=True
            else:
                i=i+1
        filein.close()

        #
        # find line in var_tdens
        #
        filepath=folder.timefun+'/time.dat'
        filein=open(str(filepath), 'r')
        input_lines = filein.readlines()
        time=float(input_lines[i+1])
        filein.close()
        click.echo(time)

    program=os.path.normpath(repo+
                             '../globals/get_data/get_data.out')

    command=(program   + ' ' +
             str(folder.result+'/tdens.dat')+ " " +
             str(time) + ' ' +
             str(file_out))
    os.system(command)

@cli.command()
@click.argument('folder_name')
# select data to compare
@click.option('flag_data', '-tdens',  flag_value='tdens',
              help='Second Option: Tdens error')
@click.option('flag_data', '-pot', flag_value='pot',
              help='Second Option: Potential error')
@click.option('flag_data', '-vel', flag_value='vel',
              help='Second Option: Velocity error')
@click.option('flag_data', '-vel_prj', flag_value='vel_prj',
              help='Second Option: Velocity error')

# select data to compare
@click.option('flag_option', '-o', flag_value='o',
              default=True,
              help='Compare w.r.t. final optimal Solution')
@click.option('flag_option', '-a', flag_value='a',
              help='Compare w.r.t. all time evolution')


def error(folder_name,flag_data,flag_option):
    """ Define new folder selecting portion of data"""
    #
    # get input files
    #
    folder=dmk.dmk_folder_structure(folder_name)
    grid_path=find_path_file(folder.abspath,'flag_grid')
    subgrid_path=find_path_file(folder.abspath,'flag_subgrid')
    
    
    if ( flag_data == 'tdens' ) :
        use_grid=grid_path
        file_exact  = str(folder.input+'/opt_tdens.dat')
        if ( flag_option =='o') :
            file_approx = str(folder.result+'/opt_tdens.dat')
        elif ( flag_option == 'a') :
            file_approx = str(folder.result+'/tdens.dat')
    
    if ( flag_data == 'vel' ) :
        use_grid=subgrid_path
        file_exact  = str(folder.input+'/opt_vel_subgrid.dat')
        if ( flag_option =='o') :
            file_approx = str(folder.result+'/opt_vel.dat')
        elif ( flag_option == 'a') :
            file_approx = str(folder.result+'/vel.dat')

    if ( flag_data == 'vel_prj' ) :
        use_grid=grid_path
        file_exact  = str(folder.input+'/opt_vel_grid.dat')
        if ( flag_option =='o') :
            file_approx = str(folder.result+'/opt_vel_prj.dat')
        elif ( flag_option == 'a') :
            file_approx = str(folder.result+'/vel_prj.dat')

    program=os.path.normpath(repo+
                             '../geometry/tools/errors_timedata.sh')
    command=(program   + ' ' +
             use_grid + ' ' +
             file_exact + ' ' +
             file_approx + ' ' +
             '1.0')
    os.system(command)

    
        
        
         
    


if __name__ == '__main__':
    cli()

