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


# Download all folders from Enrico's code
download_commands = {".": ["./init_otp_muffe.sh"]}
execute(download_commands)

download_python = {".": ["./init_py_scripts.sh"]}
execute(download_python)

build_directories = {
    "muffe_p1p0": ["./setup.sh ../Tests"],
    "muffe_sparse_optimization": ["./setup.sh simplifications"],
    "python_scripts": [
        "mv build_run_getgraph.py ../Tests",
        "mv continuous2graph.py ../Tests",
        "mv debugger.py ../Tests",
        "mv discrete2graph.py ../Tests",
        "mv filtering.py ../Tests/",
        "mv get_graph.py ../Tests/",
        "mv Getting_sources_and_sinks.py ../Tests/",
        "mv graph2plotly.py ../Tests/",
        "mv main.py ../Tests",
        "mv pre_extraction.py ../Tests",
        "mv quality_measure.py ../Tests/",
        "mv source_sink_generator.py ../Tests/",
        "mv steiner_simplification.py ../Tests/",
        "mv transport_networks.py ../Tests",
        "mv utilities.py ../Tests",
        "mv utils.py ../Tests",
        "mv terminal_computation.py ../Tests",
        "mv utilities.py ../Tests/python_script" "mv graph_cell_folder ../Tests",
        "mv par_files ../muffe_sparse_optimization/simplifications/",
        "mv images ../Tests",
        "mv PVM_data ../Tests",
    ],
}
execute(build_directories)

# remove_directories = {"python_scripts": ["rm -rf ;python_scripts"]}
# execute(remove_directories)


# Build commands
build_commands = {
    "geometry": ["./make_all.sh",],
    "globals": ["./make all.sh",],
    "linear_algebra/code": ["make clobber", "make dirs", "make",],
    "linear_algebra/Mxv": ["make clobber", "make dirs", "make",],
    "linear_algebra/test": ["make clobber", "make dirs", "make",],
    "muffe_p1p0/code": ["make clobber", "make dirs", "make",],
    "muffe_p1p0/preprocess/2d_assembly/subgrid_preprocess/code": [
        "make clobber",
        "make dirs",
        "make",
    ],
    "muffe_p1p0/python_interface": ["make clobber", "make dirs", "make",],
    "muffe_sparse_optimization/code": ["make clobber", "make dirs", "make",],
    "p1galerkin/code": ["make clobber", "make dirs", "make",],
    "p1galerkin/assembly_stiff": ["make clobber", "make dirs", "make",],
    "p1galerkin/dirac2rhs": ["make clobber", "make dirs", "make",],
    "p1galerkin/eval_divergence_constrain": ["make clobber", "make dirs", "make",],
    "p1galerkin/eval_gradient": ["make clobber", "make dirs", "make",],
    "p1galerkin/makerhs": ["make clobber", "make dirs", "make",],
    "p1galerkin/Laplace_Convergence": ["make clobber", "make dirs", "make",],
    "spectral/src": ["make clobber", "make dirs", "make",],
    "vtk": ["make clobber", "make dirs", "make",],
}
execute(build_commands)

make_command = {
    "muffe_p1p0/code": ["make",],
    "globals/axpy_timedata ": ["make",],
    "geometry/interpolate_timedata/":: ["make",],
}
execute(make_command)
