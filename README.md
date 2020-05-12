# Network extraction by routing optimization

NextRout (Network Extraction from Routing Optimization) is a tool that allows you to extract Network topologies from dynamical equations related to Routing Optimization problems. It is divided into three steps: extraction of an optimal trajectory from a dynamical system of equations in a continuous framework, pre-extraction of the network, and filtering of possible redundancies.

You can use inputs such as images, giving a cleaning that allows you to have a less redundant network. You can choose whether only have a graph from a continuous framework on routing optimization, or go beyond and apply a similar dynamics that allows you to execute a filtering without losing important information stored on your network. 

## Getting Started

### Prerequisites

The required dependencies are

---For the DMK-Continuous ---

1 - A Fortran compiler (any gfortran version>4 but version 5);
2 - Blas and Lapach libraries;
3 - Python 3;	
4 - Meshpy;
5 - PyVTK;

** You should have Python 3 installed in order to perform the following steps.

### Installation

(In this section, you have a step by step series of examples that tell you how to run a simulation).
--------------------------------------------------------------------------------------------------------------------------
Linux

** First, clone this repository by executing the following commands:

```
$ git clone https://github.com/Danielaleite/Nextrout
cd Nextrout
python setup.py
```

This is going to download the "DMK-solver" and the required python scripts, and prepare the setup for you to be able to run a simulation.
The other files inside this folder are being already executed by the following scripts (so you don't have to touch them): "init_otp_muffe.sh" takes all the repositories at https://gitlab.com/opt-muffe/otp_utilities, while "init_py_scripts.sh" takes the required scripts at https://github.com/Danielaleite/python_scripts.  


The execution of such scripts takes around 8 minutes to finish (so dont't panic if it's taking a while to be done). In the end, you should have the following folders:

- globals
-linear algebra
- geometry
- muffe_p1p0
- muffe_sparse_optimization
- p1galerkin
- python_scripts
- spectral
- Tests
- vtk
------------------------------------------------------------

## Running the tests

To be able to run an experiment, go to the folder "Tests" and run the script "build_run_get_graph.py".

```
cd network_extraction
python network_extraction.py
```

The following questions are going to appear:

* What is your flag?
(add answer to this question)
* What is the number of divisions for the Mesh?
* What is the value for beta?
* What is(are) the flag(s) for source?
* What is(are) the flag(s) for sink?
* What is the val for the discrete beta?: 
* continuous DMK solver?: 
* graph extraction?
* graph simplification, reduction?

If everything goes well, you are going to have the results located inside the folder "muffe_sparse_optimization/simplifications/runs/"
 

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






""""""""""""""""""""""""OLD""""""""""""""""""""""""""""""""""""

# Network extraction by routing optimization

NextRout (Network Extraction from Routing Optimization) is a tool that allows you to extract Networks from different inputs, such as images, giving a cleaning that allows you to have a less redundant network. You can choose whether only have a graph from a continuous framework on routing optimization, or go beyond and apply a similar dynamics that allows you to execute a filtering without losing important information stored on your network. 

## Getting Started

### Prerequisites

The required dependencies are

---For the DMK-Continuous ---

1 - A Fortran compiler (any gfortran version>4 but version 5);
2 - Blas and Lapach libraries;
3 - Python 3;
4 - Meshpy;
5 - PyVTK;

** You should have Python 3 installed in order to perform the following steps.

### Installation

(In this section, you have a step by step series of examples that tell you how to run a simulation).
--------------------------------------------------------------------------------------------------------------------------
Linux

** First, clone this repository by executing the following commands:

```
$ git clone https://github.com/Danielaleite/Netwok_extraction_from_continuous
cd Netwok_extraction_from_continuous
python setup.py
```

This is going to download the "DMK-solver" and python scripts, and prepare the setup for you to be able to run a simulation.
The other files inside this folder are being already executed by the following scripts (so you don't have to touch them): "init_otp_muffe.sh" takes all the repositories at https://gitlab.com/opt-muffe/otp_utilities, while "init_py_scripts.sh" takes the required scripts at https://github.com/Danielaleite/python_scripts.  


The execution of such scripts takes around 8 minutes to finish (so dont't panic if it's taking a while to be done). In the end, you should have the following folders:

-------- Give a description for each folder ---------------

- globals
-linear algebra
- geometry
- muffe_p1p0
- muffe_sparse_optimization
- p1galerkin
- python_scripts
- spectral
- Tests
- vtk
------------------------------------------------------------

## Running the tests

To be able to run an experiment, go to the folder "Tests" and run the script "build_run_get_graph.py".

```
cd Tests
python build_run_get_graph.py
```

The following questions are going to appear:

* What is your flag?
(add answer to this question)
* What is the number of divisions for the Mesh?
* What is the value for beta?
* What is(are) the flag(s) for source?
* What is(are) the flag(s) for sink?
* What is the val for the discrete beta?: 
* continuous DMK solver?: 
* graph extraction?
* graph simplification, reduction?

If everything goes well, you are going to have the results located inside the folder "muffe_sparse_optimization/simplifications/runs/"
 

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


"""