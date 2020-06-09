# Network extraction by routing optimization

**NextRout** - Network Extraction from Routing Optimization - is a tool created to extract Network topologies from dynamical equations related to routing optimization problems. It is divided into three steps: extraction of an optimal trajectory from a dynamical system of equations in a continuous framework, pre-extraction of the network, and filtering of possible redundancies.

You can use in input pre-build functions, graphs or images.  The filtering allows you to have a less redundant network. You can choose whether only have a graph from a continuous framework on routing optimization, or go beyond and apply a similar dynamics that allows to execute a filtering without losing important information stored on the network. 

## Prerequisites

The required dependencies are

* For the DMK-Continuous 

	* A Fortran compiler (any gfortran version>4 but version 5);
 	* Blas and Lapach libraries;
 	* Python 3;	
 	* Meshpy;
 	* PyVTK;

* For the DMK-discrete
	* **Python 3** is enough.

## Installation

**Linux/OSX**

1. Clone this repository:

```
git clone https://github.com/Danielaleite/Nextrout
cd Nextrout
python setup.py
```

This is going to download the "DMK-solver" and all required python modules. The execution of "setup.py" takes around 9 minutes to be finished (so dont't panic if you notice it's taking a while). In the end, you should have the following folders:

* otp_utilities
* python_scripts
* results

Inside otp_utilities you will see a bash script "init_otp_muffe.sh". This is a piece previously executed in "setup.py" that takes all the repositories at https://gitlab.com/opt-muffe/otp_utilities, so don't worry about it. After executing setup, you will be able to see three folders,

* **otp_utilities**
* **python_scripts**
* **results**

You will find in each folder another README.me with a more detailed description to the subfolders and main scripts. 


## How to perform a simulation?

Inside "python_scripts" look for "network_extraction.py" ,

```
cd python_scripts
python network_extraction.py
```

For the input parameters, you should reply the following questions:

1. Do you have an input.ctrl file? [y/N]:
If your answer is yes, skip to the question 8. If not, 
2.  What is the input flag?
This is where you give the initial transport density. There are three options: predefined funcions for DMK-continuous (tdens0), images or your own graph; for more details, see https://github.com/Danielaleite/Nextrout/blob/master/python_scripts/default_inputs.md.
3. What is the number of divisions of the mesh?
4. Which beta should be used in the DMK-solver?
5. What is are) the flag(s) for source function?
6. What is(are) the flag(s) for sink function?
7. Which beta should be used in the filtering?
8. Should the DMK-solver be executed?
9. Do you want the graph pre-extraction to be done?
10. Do you need the graph filtering/reduction to be done?


After this, depending of your parameters, it might from a few seconds to a few minutes untill the simulation is performed. Th results are all stored inside the folder "results";

## Contributing

Nextrout is being developed continuously. We appreciate any help or suggestions. Do not hesitate to contact us in case you want to contribute or have any questions.

## Versioning


## Authors

* Daniela Leite, Diego Theuerkauf 

See also the list of [contributors](https://github.com/Danielaleite/Nextrout/graphs/contributors) who participated in this project.