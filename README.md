# Network extraction by routing optimization

**NextRout** - Network Extraction from Routing Optimization - is a tool created to extract Network topologies from dynamical equations related to routing optimization problems. It is divided into three steps: extraction of an optimal trajectory from a dynamical system of equations in a continuous framework; pre-extraction of the network; and filtering of possible redundancies.

You can use in input pre-build functions, graphs or images.  The filtering allows you to have a less redundant network. You can choose whether only have a graph from a continuous framework on routing optimization, or go beyond and apply a similar dynamics that allows to execute a filtering without losing important information stored on the network. 

NextRout is based on the theory described in this paper:

- _Network extraction by routing optimization_ [published here](https://www.nature.com/articles/s41598-020-77064-4), D. Baptista, D. Leite, E. Facca, M. Putti and C. D. Bacco, arxiv 2005.02805, 2020.   
If you use this software please cite this work.

## Prerequisites

The required dependencies are

* For the DMK-Continuous 

	* A Fortran compiler (any gfortran version>4 but version 5);
 	* Blas and Lapack libraries;
 	* Python 3;	
 	* Meshpy;
 	* PyVTK;

* For the DMK-discrete
	* **Python 3** is enough.

## Installation

**Linux/OSX**

Just clone this repository:

```
git clone https://github.com/Danielaleite/Nextrout
cd Nextrout
python setup.py
```

This is going to download the "DMK-solver" and all required python dependencies. The execution of "setup.py" takes around 9 minutes to be finished (so dont't panic if you notice it's taking a while). In the end, you should have the following folders:

* **otp_utilities**
* **python_scripts**
* **results**

Inside **otp_utilities** there is a bash script named "init_otp_muffe.sh". This is a piece previously executed in "setup.py" that takes all the repositories at https://gitlab.com/opt-muffe/otp_utilities (within this link you can find a more detailed description to the _DMK-solver_). In **python_scripts**, as the name suggests, you will be able to find the main python scripts and subfolders where the entire procedure takes place. The folder **results** stores the final graphs and solutions from *DMK-continuous*. You can find in each folder another README.md with a more detailed description and instructions for each subfolder/script. 


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


After this, depending of your input parameters, it might take from a few seconds to minutes untill the entire procedure is completely done. To see the results, just do

```
cd results
```

## Contributing

Nextrout is being developed continuously. We appreciate any help or suggestions. Do not hesitate to contact us in case you want to contribute or have any questions.


## Authors

* Daniela Leite, Diego Theuerkauf 

See also the list of [contributors](https://github.com/Danielaleite/Nextrout/graphs/contributors) who participated in this project.

## License

Copyright (c) 2020 Daniela Leite, Diego Theuerkauf and Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

