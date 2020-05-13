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
"""
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
"""
build_directories = {
    "otp_utilities/muffe_p1p0": ["./setup.sh ../../python_scripts"],
    "otp_utilities/muffe_sparse_optimization/": ["./setup.sh simplifications"],
    "python_scripts": ["cp utilities.py ./python_script",],
}
"""
        "mv continuous2graph.py ../python_scripts",
        "mv debugger.py ../python_scripts",
        "mv discrete2graph.py ../python_scripts",
        "mv filtering.py ../python_scripts",
        "mv get_graph.py ../python_scripts",
        "mv Getting_sources_and_sinks.py ../python_scripts",
        "mv graph2plotly.py ../python_scripts",
        "mv main.py ../python_scripts",
        "mv pre_extraction.py ../python_scripts",
        "mv quality_measure.py ../python_scripts",
        "mv source_sink_generator.py ../python_scripts",
        "mv steiner_simplification.py ../python_scripts",
        "mv transport_networks.py ../python_scripts",
        "mv utilities.py ../python_scripts",
        "mv utils.py ../python_scripts",
        "mv terminal_computation.py ../python_scripts",
        "mv utilities.py ../python_scripts/python_script",
        "mv graph_cell_folder ../python_scripts/",
        "mv par_files ../python_scripts/simplifications/",
        "mv images ../python_scripts",
        "mv PVM_data ../python_scripts",
    ],
"""


execute(build_directories)

"""
make_command = {
    "otp_utilities/muffe_p1p0/code": ["make",],
    "otp_utilities/globals/axpy_timedata": ["make",],
    "otp_utilities/geometry/interpolate_timedata/": ["make",],
}
execute(make_command)
"""
