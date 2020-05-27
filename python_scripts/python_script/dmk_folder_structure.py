import os
import utilities as util


class dmk_folder_structure:
    #
    # set default structure
    #
    def __init__(self,folder_path):
        #
        # set name local path and abssolute folders
        #
        self.path         = os.path.normpath(folder_path)
        self.name         = os.path.basename(folder_path)
        self.directory    = os.path.abspath(os.path.dirname(folder_path))
        self.abspath      = os.path.abspath(folder_path)
        
        #
        # set subfolders names
        #
        self.input        = self.abspath   + '/input'
        self.output       = self.abspath   + '/output'
        self.input_vtk    = self.input  + '/vtk'
        self.output_vtk   = self.output + '/vtk'
        self.result       = self.output + '/result'
        self.timefun      = self.output + '/timefun'
        self.linsys       = self.output + '/linsys'
        self.figures      = self.output + '/figures'
        
    
    def mkdirs(self):
        # 
        # create subfolders
        #
        if( not os.path.exists(self.abspath) ) :
            os.mkdir(self.abspath)
            os.mkdir(str(self.input))
            os.mkdir(str(self.input_vtk))
            os.mkdir(str(self.output))
            os.mkdir(str(self.result))
            os.mkdir(str(self.timefun))
            os.mkdir(str(self.linsys))
            os.mkdir(str(self.output_vtk))
    
