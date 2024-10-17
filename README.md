<p align="center">
  <img src="https://img.shields.io/badge/python-3.6-blue.svg" alt="Python 3.6">
</p>

---
# Network extraction by routing optimization

**NextRout** - Network Extraction from Routing Optimization - is a tool created to extract Network topologies from dynamical equations related to routing optimization problems. It is divided into three steps: extraction of an optimal trajectory from a dynamical system of equations in a continuous framework; pre-extraction of the network; and filtering of possible redundancies.

You can use in input pre-build functions, graphs or images.  The filtering allows you to have a less redundant network. You can choose whether only have a graph from a continuous framework on routing optimization, or go beyond and apply a similar dynamics that allows to execute a filtering without losing important information stored on the network. 

NextRout is based on the theory described in this paper:

- [_Network extraction by routing optimization_](https://www.nature.com/articles/s41598-020-77064-4), D. Baptista, D. Leite, E. Facca, M. Putti and C. D. Bacco, Nature Scientific Reports, 10, 20806, 2020.   
If you use this software please cite this work.

## Prerequisites

The required dependencies are

	* A Fortran compiler (any gfortran version>4 but version 5);
 	* Blas and Lapack libraries;
 	* Python 3;	
 	* Meshpy;
    * Click;
    * Numpy v<1.19;
    * f90wrap; 

## Installation

**Linux/OSX**

Once all the dependencies are installed, you just need to clone this repository:

```
git clone https://github.com/Danielaleite/Nextrout
cd Nextrout
python setup.py
```

This is going to download the "DMK-solver" and all required python dependencies.  

The script `setup.py` is using the default Fortran compile `/usr/bin/gfortran`. 
If you want to use another path for this compiler, then change the line in side `setup.py`:
```bash
os.system('cd dmk_utilities/dmk_solver && mkdir build/ &&  cd build &&  cmake ../ &&  make')
```
to
```
os.system('cd dmk_utilities/dmk_solver && mkdir build/ &&  cd build && cmake -DMy_Fortran_Compiler=your_fortran_compiler_path ../ &&  make')
```
where `your_fortran_compiler_path` is your custom path, e.g. `/usr/local/bin/gfortran`.

The execution of `setup.py` takes a few minutes to be finished (so dont't panic if you notice it's taking a while). In the end, you should have the following folders:

* **otp_utilities**
* **nextrout_core**

Inside **otp_utilities** there are all the files related to the DMK solver. These are needed to execute **NextRout**. In **nextrout_core**, as the name suggests, you will be able to find the main python scripts and subfolders where the entire procedure takes place. In case of errors during installation, please visit [_DMK solver_](https://gitlab.com/enrico_facca/dmk_solver), section **Troubleshooting**. 


## How to perform a simulation?

Inside "nextrout_core/" look for "test.py" ,

```
cd nextrout_core
python test.py
```

Results will be stored in the folder ```outputs.```

## Contributing

Nextrout is being developed continuously. We appreciate any help or suggestions. Do not hesitate to contact us in case you want to contribute or have any questions.


## Authors

* Daniela Leite, Diego Theuerkauf, Caterina De Bacco, Enrico Facca, Marco Putti

See also the list of [contributors](https://github.com/Danielaleite/Nextrout/graphs/contributors) who participated in this project.

## License

Copyright (c) 2020 Daniela Leite, Diego Theuerkauf and Caterina De Bacco

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

