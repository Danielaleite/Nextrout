# Network extraction by routing optimization

NextRout (Network Extraction from Routing Optimization) is a tool that allows you to extract Network topologies from dynamical equations related to Routing Optimization problems. It is divided into three steps: extraction of an optimal trajectory from a dynamical system of equations in a continuous framework, pre-extraction of the network, and filtering of possible redundancies.

You can use in input pre-build functions, graphs or images.  The filtering allows you to have a less redundant network. You can choose whether only have a graph from a continuous framework on routing optimization, or go beyond and apply a similar dynamics that allows to execute a filtering without losing important information stored on the network. 

## Getting Started

### Prerequisites

The required dependencies are

---For the DMK-Continuous ---

1 - A Fortran compiler (any gfortran version>4 but version 5);
2 - Blas and Lapach libraries;
3 - Python 3;	
4 - Meshpy;
5 - PyVTK;

** You should have Python 3 installed in order to perform the next steps.

### Installation

Linux

** First, clone this repository by executing the following commands:

```
git clone https://github.com/Danielaleite/Nextrout
cd Nextrout
python setup.py
```

This is going to download the "DMK-solver" and the required python modules and prepare the setup for running a simulation. The execution of "setup.py" takes around 8 minutes to be finished (so dont't panic if you notice it's taking a while). In the end, you should have the following folders:

* otp_utilities
* python_scripts
* results

Inside otp_utilities you will find "init_otp_muffe.sh". This is a piece previously executed in "setup.py" that takes all the repositories at https://gitlab.com/opt-muffe/otp_utilities, so you shouldn't worry about it. After executing setup, you will be able to see the following respositories:

* globals
*linear algebra
* geometry
* muffe_p1p0
* muffe_sparse_optimization
* p1galerkin
* python_scripts
* spectral
* Tests
- vtk

In "python_scripts" you should have the following content:

* PVM_data
* runs
* network_impainting
* network_analysis
* bash_scripts
* python_script
* jobs
* images
* graph_cell_folder


* build_run_getgraph.py
* continuous2graph.py
* debugger.py
* discrete2graph.py
* dmk_folder.py
* filtering.py
* fort.555
* get_graph.py
* Getting_sources_and_sinks.py
* graph2plotly.py
* inputs.ctrl
* Instructions.md
* location_repository.txt
* main.py
* muffa.ctrl
* network_extraction.py
* new_muffa.ctrl
* outputs_dmk_c.txt
* outputs_gg.txt
* outputs_vtk.txt
* pre_extraction.py
* quality_measure.py
* README.md
* setup.py
* source_sink_generator.py
* steiner_simplification.py
* terminal_computation.py
* transport_networks.py
* utilities.py
* utils.py

Inside "results" you will find two subfolders that are going to be filled with the results of the experiments. 

* discrete
* continuous
------------------------------------------------------------

## How to perform a simulation?

Inside "python_scripts" you should run "network_extraction.py". So,

```
cd python_scripts
python network_extraction.py
```

You will see the following questions:

* What is the tdens0 flag?
This is where you give the initial transport density. You have three options: predefined funcions for DMK-continuous, images or your own graph;
* What is the number of divisions of the mesh?

* What is the value for beta to be used in the DMK-solver?
* What is are) the flag(s) for source function?
* What is(are) the flag(s) for sink function?
* What is the value for beta to be used in the filtering part?
* Should the DMK-solver step be executed?
* Should the graph pre-extraction be done?
* Should the graph filtering/reduction be done? 


 

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```


## Built With


## Contributing

## Versioning


## Authors

* **Daniela Leite** **Diego Theuerkauf** 

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments
