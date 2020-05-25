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

click.echo("Define the parameters for the network extraction routine.")


@click.command()
@click.option(
    "--flag_list",
    prompt="What is the tdens0 flag?",
    help="This is where you should put the name of the tdens0 flag.",
)
@click.option(
    "--ndiv_list",
    prompt="What is the number of divisions of the mesh?",
    help="This is the number of divisions for your triangular mesh.",
)
@click.option(
    "--beta_list",
    prompt="What is the value for beta to be used in the DMK-solver?",
    help="This is the value of the Beta exponent.",
)
@click.option(
    "--source_list",
    prompt="What is(are) the flag(s) for source function?",
    help="Source information.",
)
@click.option(
    "--sink_list",
    prompt="What is(are) the flag(s) for sink function?",
    help="Sink information.",
)
@click.option(
    "--bd_list",
    prompt="What is the value for beta to be used in the filtering part?",
    help="Beta discrete information.",
)
@click.option(
    "--dmk_input",
    prompt="Should the DMK-solver step be executed?",
    help="Execute the DMK-solver.",
)
@click.option(
    "--ge_input",
    prompt="Should the graph pre-extraction be done?",
    help="Graph pre-extraction from the DMK solutions.",
)
@click.option(
    "--gs_input",
    prompt="Should the graph filtering/reduction be done?",
    help="Filtering the graph using discrete DMK solver. Graph reduction.",
)
def network_extraction(
    flag_list,
    beta_list,
    ndiv_list,
    source_list,
    sink_list,
    bd_list,
    dmk_input,
    ge_input,
    gs_input,
):
    """
    This script runs the network extraction procedure: DMK-solver, graph pre-extraction and graph-filtering. Input
    values must be given as "flag1,flag2,flag3,..." if many flags are used.

    :param flag_list: flags used for the initial transport density.
    :param beta_list: values for exponent beta of DMK-solver (1st step).
    :param ndiv_list: number of divisions of the x-axis used to the creation of the mesh.
    :param source_list: flags for source function.
    :param sink_list: flags for sink function. If "=", then it takes the same value as source_list.
    :param bd_list: values for exponent beta of filtering (3rd step).
    :param dmk_input: "yes" if DMK-solver to be used.
    :param ge_input: "yes" if graph pre-extraction to be done.
    :param gs_input: "yes,yes" if graph filtering *and* graph reduction to be done.
    :return:
        outputs stored in ./runs/folder_name

    """

    if sink_list == "=":
        sink_list = source_list
    source_list = source_list.split(",")
    sink_list = sink_list.split(",")
    source_sink_list = [[source_list[i], sink_list[i]] for i in range(len(source_list))]
    print(source_sink_list)
    for flag in flag_list.split(","):
        for beta in beta_list.split(","):
            for ndiv in ndiv_list.split(","):
                for source_sink in source_sink_list:
                    for beta_discr in bd_list.split(","):
                        folder_name = (
                            "%s" % flag
                            + "_"
                            + "b"
                            + str(int(10 * (float(beta))))
                            + "_"
                            + "%s" % ndiv
                            + "dv"
                            + "_sf"
                            + str(int(10 * (float(beta_discr))))
                            + "_"
                            + "%s" % source_sink[0]
                            + "%s" % source_sink[1]
                        )

                        if dmk_input == "yes":
                            click.echo("Your .ctrl file is ready! ")
                            # click.echo('Your flag is %s!' %flag)
                            file = open("inputs.ctrl", "w+")

                            file.write(
                                "# DOMAIN ####################################################### "
                                + "\n"
                            )
                            file.write("2d ! flag_domain (2d, 3d, surface)" + "\n")

                            file.write(
                                "# MESH #########################################################"
                                + "\n"
                            )
                            file.write("rect_cnst path/grid.dat ! flag_grid" + "\n")
                            file.write("%s" % ndiv + " ! ndiv" + "\n")
                            file.write(
                                "0 " + " ! nref" + "\n"
                            )  # --------------------------- !!!!!!!!!!

                            file.write(
                                "# PREBUILD EXAMPLES ############################################"
                                + "\n"
                            )
                            file.write("%s" % beta + " extra_path ! flag_pflux " + "\n")
                            file.write("1.0 extra_path ! flag_pmass (gamma)" + "\n")
                            file.write("1.0 extra_path ! flag_decay" + "\n")
                            file.write("1.0 extra ! decay0" + "\n")
                            file.write(
                                "rect_cnst  path/frog_source.dat ! flag_source" + "\n"
                            )  # %s' %source_sink[0]+
                            file.write(
                                "rect_cnst  path/frog_sink2.dat ! flag_sink" + "\n"
                            )  # %s' %source_sink[1]+
                            file.write("0 ! flag_normalize" + "\n")
                            # file.write('%s' %flag + '!flag_tdens0 ' + "\n" )
                            file.write(
                                "%s" % flag
                                + " "
                                + "path/tdens0.dat ! flag_tdens0 "
                                + "\n"
                            )
                            file.write("1.0 extra_path ! flag_kappa" + "\n")

                            file.write(
                                "## TIME RANGE ###################################################"
                                + "\n"
                            )
                            file.write("0.0 ! tzero" + "\n")
                            file.write("5.0e2 ! tmax (stop while loop)")
                            # list = file.readlines()
                            file.close()

                            # folder_name =  "%s" %flag + "_" + "b"  +str(int(10*(float(beta)))) + '_' + "%s" %ndiv + "dv" + '_sf'+str(int(10*(float(beta_discr))))+ '_'+ "%s" %source_sink[0]+ "%s" %source_sink[1]

                            def init(folder_name):
                                new_dir = "./runs/" + folder_name
                                try:
                                    os.mkdir(new_dir)
                                except OSError:
                                    print(
                                        "Creation of the directory %s failed." % new_dir
                                    )

                            try:
                                shutil.rmtree("./runs/" + folder_name)
                            except OSError:
                                pass
                            command = (
                                "./dmk_folder.py assembly "
                                + "./runs/"
                                + folder_name
                                + " inputs.ctrl "
                            )
                            os.system(command)

                            source_sink_generator(
                                "./runs/" + folder_name,
                                ndiv,
                                source_sink[0],
                                source_sink[1],
                            )
                            if (
                                source_sink[0] != "rect_cnst"
                                and source_sink[1] != "rect_cnst"
                            ):
                                source_sink_preprocess("./runs/" + folder_name)

                            command = (
                                "./dmk_folder.py run "
                                + "./runs/"
                                + folder_name
                                + " new_muffa.ctrl > outputs_dmk_c.txt"
                            )
                            os.system(command)

                            command = (
                                "./dmk_folder.py get-graph "
                                + "./runs/"
                                + folder_name
                                + " 0.1 "
                                + " 100000000000  > outputs_gg.txt"
                            )
                            os.system(command)

                            command = (
                                "./dmk_folder.py vtk -a -tdens "
                                + "./runs/"
                                + folder_name
                                + " > outputs_vtk.txt"
                            )
                            os.system(command)

                        else:
                            print("Skipping DMK-solver part.")

                        ############## GRAPH PRE-EXTRACTION######################

                        errors = []
                        for funct in ["tdens"]:
                            print("Now we are running: ", funct)

                            for graph_type in ["1"]:  # ,'2','3']:
                                print("Now we are running graph: ", graph_type)
                                for weighting_method_graph in ["ER"]:  # , 'AVG']:
                                    # print('Now we are running weighting method graph: ', weighting_method_graph)
                                    if funct == "tdens":
                                        if 0 < float(beta) <= 1:
                                            t_list = [0.1]  # ,.65]
                                        else:
                                            t_list = [0.001]  # <------------ !!
                                    else:  # flux
                                        if 0 < float(beta) <= 1:
                                            t_list = [0.1]  # ,.65]
                                        else:
                                            t_list = [0.001]

                                    for threshold in t_list:
                                        subfolder = "./runs/" + folder_name
                                        t = float(threshold)

                                        new_dir = subfolder + "/" + funct

                                        try:
                                            os.mkdir(new_dir)
                                        except OSError:
                                            print(
                                                "Creation of the directory %s failed."
                                                % new_dir
                                            )

                                        if ge_input == "yes":
                                            print(
                                                "Flag and tolerance:",
                                                subfolder,
                                                t,
                                                graph_type,
                                                funct,
                                            )
                                            G = graph_extraction_from_dat_files(
                                                subfolder,
                                                t,
                                                graph_type,
                                                funct,
                                                weighting_method_graph,
                                                source_sink[0],
                                                source_sink[1],
                                            )
                                            ### writing parameters ###
                                            file = open(
                                                subfolder + "/parameters.ctrl", "w+"
                                            )
                                            file.write(
                                                "--- dmk configuration ---" + "\n"
                                            )
                                            file.write("tdens0: " + str(flag) + "\n")
                                            file.write("ndiv: " + str(ndiv) + "\n")
                                            file.write(
                                                "source_flag: "
                                                + str(source_sink[0])
                                                + "\n"
                                            )
                                            file.write(
                                                "sink_flag: "
                                                + str(source_sink[1])
                                                + "\n"
                                            )
                                            file.write("beta: " + str(beta) + "\n")
                                            file.write(
                                                "--- graph pre-extraction ---" + "\n"
                                            )
                                            file.write(
                                                "delta: " + str(threshold) + "\n"
                                            )
                                            file.write(
                                                "graph_type: " + str(graph_type) + "\n"
                                            )
                                            file.write("funct: " + str(funct) + "\n")
                                            file.write(
                                                "weighting_method_graph: "
                                                + str(weighting_method_graph)
                                                + "\n"
                                            )

                                        else:
                                            print("Skipping graph-extraction part.")
                                        print(gs_input.split(","))
                                        if gs_input.split(",")[0] == "yes":
                                            for minimum in [0.001]:  # [0.001]:
                                                for weighting_method_simplification in [
                                                    "ER"
                                                ]:  # IBP, BPW
                                                    for btns_factor in [
                                                        (-1, -1)
                                                    ]:  # 0.010]
                                                        print(
                                                            "=======================================================",
                                                            graph_type,
                                                            funct,
                                                            weighting_method_graph,
                                                            weighting_method_simplification,
                                                        )
                                                        min_ = float(minimum)
                                                        btns_factor_source = float(
                                                            btns_factor[0]
                                                        )
                                                        btns_factor_sink = float(
                                                            btns_factor[1]
                                                        )
                                                        BP_weights = "BPtdens"
                                                        reduction_flag = gs_input.split(
                                                            ","
                                                        )[
                                                            1
                                                        ]  # write 'yes' to get the 2nd simpl
                                                        # print('>>>>>____________________Computing bp simplification for',subfolder, funct, t, graph_type, min_,weighting_method_graph, weighting_method_simplification)
                                                        i = 0
                                                        beta_discr = float(beta_discr)
                                                        i += 1

                                                        graph_filtering_from_dat_files(
                                                            subfolder,
                                                            t,
                                                            graph_type,
                                                            beta_discr,
                                                            funct,
                                                            min_,
                                                            btns_factor_source,
                                                            btns_factor_sink,
                                                            weighting_method_graph,
                                                            weighting_method_simplification,
                                                            source_sink[0],
                                                            source_sink[1],
                                                            BP_weights,
                                                            reduction_flag,
                                                        )

                                                        file.write(
                                                            "--- filtering ---" + "\n"
                                                        )
                                                        file.write(
                                                            "beta: "
                                                            + str(beta_discr)
                                                            + "\n"
                                                        )
                                                        file.write(
                                                            "delta: " + str(min_) + "\n"
                                                        )
                                                        file.write(
                                                            "btns_factor (s+/s-): "
                                                            + str(btns_factor)
                                                            + "\n"
                                                        )
                                                        file.write(
                                                            "edge weights built from: "
                                                            + str(BP_weights)
                                                            + "\n"
                                                        )
                                                        file.write(
                                                            "degree-2 reduction: "
                                                            + str(reduction_flag)
                                                            + "\n"
                                                        )
                        dest = "../results/continuous"
                        shutil.move("runs/" + folder_name, dest)

                        dest2 = "../results/discrete"
                        shutil.move(
                            "../otp_utilities/muffe_sparse_optimization/simplifications/runs/"
                            + folder_name,
                            dest2,
                        )


##############

if __name__ == "__main__":
    network_extraction()
