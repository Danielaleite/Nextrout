#!usr/bin/env python3

import os
import os.path
import contextlib


@contextlib.contextmanager
def change_temporarily_into(dir):
    """doing cd into dir and go back to the original directory afterwards."""
    current_dir = os.getcwd()

    try:
        os.chdir(dir)
        yield
    except Exception as e:
        if isinstance(e, FileNotFoundError):
            print("I am in {}".format(os.getcwd()))
            print("You are trying to cd into {}".format(dir))
            print("Available directories:{}".format(os.listdir(os.getcwd())))
        raise e
    finally:
        os.chdir(current_dir)


def execute(commands_in_dir):
    """Execute a list of commands in a directory.
     
     Keys are directories, values are lists of commands.
    """

    def execute_in_dir(dir, commands):
        with change_temporarily_into(dir):
            for cmd in commands:
                os.system(cmd)

    for dir, commands in list(commands_in_dir.items()):
        execute_in_dir(dir, commands)


# Download dependencies for the first step and backend for filtering
download_commands = {"./otp_utilities": ["./init_otp_muffe.sh"]}
execute(download_commands)

# Build commands
build_commands = {
    "otp_utilities/geometry": ["./make_all.sh",],
    "otp_utilities/globals": ["./make all.sh",],
    "otp_utilities/linear_algebra/code": ["make clobber", "make dirs", "make",],
    "otp_utilities/linear_algebra/Mxv": ["make clobber", "make dirs", "make",],
    "otp_utilities/linear_algebra/test": ["make clobber", "make dirs", "make",],
    "otp_utilities/muffe_p1p0/code": ["make clobber", "make dirs", "make",],
    "otp_utilities/muffe_p1p0/preprocess/2d_assembly/subgrid_preprocess/code": [
        "make clobber",
        "make dirs",
        "make",
    ],
    "otp_utilities/muffe_p1p0/python_interface": ["make clobber", "make dirs", "make",],
    "otp_utilities/muffe_sparse_optimization/code": [
        "make clobber",
        "make dirs",
        "make",
    ],
    "otp_utilities/p1galerkin/code": ["make clobber", "make dirs", "make",],
    "otp_utilities/p1galerkin/assembly_stiff": ["make clobber", "make dirs", "make",],
    "otp_utilities/p1galerkin/dirac2rhs": ["make clobber", "make dirs", "make",],
    "otp_utilities/p1galerkin/eval_divergence_constrain": [
        "make clobber",
        "make dirs",
        "make",
    ],
    "otp_utilities/p1galerkin/eval_gradient": ["make clobber", "make dirs", "make",],
    "otp_utilities/p1galerkin/makerhs": ["make clobber", "make dirs", "make",],
    "otp_utilities/p1galerkin/Laplace_Convergence": [
        "make clobber",
        "make dirs",
        "make",
    ],
    "otp_utilities/spectral/src": ["make clobber", "make dirs", "make",],
    "otp_utilities/vtk": ["make clobber", "make dirs", "make",],
}

execute(build_commands)

build_directories = {
    "otp_utilities/muffe_p1p0": ["./setup.sh ../../network_extraction"],
    "otp_utilities/muffe_sparse_optimization/": [
        "./setup.sh ../../network_extraction/simplifications"
    ],
    "python_scripts": [
        "cp network_extraction.py ../network_extraction",
        "cp continuous2graph.py ../network_extraction",
        "cp debugger.py ../network_extraction",
        "cp discrete2graph.py ../network_extraction",
        "cp filtering.py ../network_extraction",
        "cp get_graph.py ../network_extraction",
        "cp Getting_sources_and_sinks.py ../network_extraction",
        "cp graph2plotly.py ../network_extraction",
        "cp main.py ../network_extraction",
        "cp pre_extraction.py ../network_extraction",
        "cp quality_measure.py ../network_extraction",
        "cp source_sink_generator.py ../network_extraction",
        "cp steiner_simplification.py ../network_extraction",
        "cp transport_networks.py ../network_extraction",
        "cp utilities.py ../network_extraction",
        "cp utils.py ../network_extraction",
        "cp terminal_computation.py ../network_extraction",
        "cp utilities.py ../network_extraction/python_script",
        "cp -a graph_cell_folder ../network_extraction/",
        "cp -a par_files ../network_extraction/simplifications/",
        "cp -a images ../network_extraction",
        "cp -a PVM_data ../network_extraction",
    ],
}

execute(build_directories)

# remove_directories = {"python_scripts": ["rm -rf ;python_scripts"]}
# execute(remove_directories)

"""
make_command = {
    "muffe_p1p0/code": ["make",],
    "globals/axpy_timedata": ["make",],
    "geometry/interpolate_timedata/": ["make",],
}
execute(make_command)
"""
