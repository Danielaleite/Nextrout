import os
from filtering_triangles import *


def get_many_graphs(folder_list, threshold):
    for folder_name in folder_list:
        folder_name = "runs/" + folder_name
        opt_tdens = extracting_weights(folder_name, "output/result/opt_tdens.dat")
        _, weights_tdens = weight2dict(opt_tdens, "output/result/opt_tdens.dat")
        maxW = max(weights_tdens)
        command = (
            "./dmk_folder.py get-graph "
            + folder_name
            + " "
            + str(maxW * threshold)
            + " 1000"
        )
        os.system(command)


folder_list = [
    "1_b11_6dv_4rcl4rcl",
    "1_b11_6dv_4rcm4rcm",
    "1_b11_6dv_5rcl5rcl",
    "1_b11_6dv_5rch5rch",
    "1_b11_6dv_5rcm5rcm",
    "1_b11_6dv_2rcl2rcl",
    "1_b11_6dv_2rcs2rcs",
    "1_b11_6dv_3rch3rch",
    "1_b11_6dv_3rcl3rcl",
    "1_b11_6dv_3rc3rc",
]
folder_list = ["1_b11_6dv_5rch5rch_1"]
get_many_graphs(folder_list, 0.01)
