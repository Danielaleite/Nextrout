import click
import os
import shutil
import distutils
import numpy as np
import shutil
import distutils
from threading import Thread
from filtering_triangles import *
from BP_simplification import *
from source_sink_generator import *
#import BP_simplification
import re
import networkx as nx
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt



import os, os.path
import shutil
import glob
import subprocess
from shutil import copytree, ignore_patterns


click.echo('Define the parameters for the continuous part here!')
@click.command()
@click.option('--source_list', prompt='What is your source?', help='This is where you should write the name of the source flag.')
@click.option('--sink_list', prompt='What is your sink?', help='This is where you should write the name of the sink flag.')
@click.option('--transport_cont_list', prompt='What is your sink?', help='This is where you should write the name of the sink flag.')
@click.option('--transport_discrete_list', prompt='What is your sink?', help='This is where you should write the name of the sink flag.')

def transport(flag_list,beta_list,ndiv_list,source_list,sink_list):
