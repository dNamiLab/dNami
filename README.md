# dNami
**dNami** [di:nɑ:mi:] is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes on structured meshes. 
dNami uses MPI and cache blocking techniques to speed up stencil-based operations.

<p align="center">
  <img src="./doc/usage/img/earth.gif" alt="earth.gif" />
</p>

<p align="center">
	<i> Animation of the atmospheric and water-height disturbance due to the Tonga volcano explosion in January 2022 computed with dNami. Check out a higher resolution version on YouTube: https://www.youtube.com/watch?v=sIivtcx1Al0</i>
</p>

## Motivation and scope

The time evolution of a variety of physical and biological processes may be described by systems of balance laws, which, if given appropriate initial and boundary conditions, dictate the future states of the systems. For instance, systems of balance laws invoking mass, momentum and energy have been incredibly successful at providing meaningful insights to the future states of realistic systems in physics (e.g. fluid dynamics). Yet, experimenting numerically with such systems still requires much implementation time. dNami  was created so that more research time is spent exploring the dynamical properties of the system of balance laws of interest to the user, and less time is wasted on its numerical implementation across the whole computational spectrum, from the initial small-scale exploratory work on a workstation to the final large-scale computations on national clusters. Thus, dNami is a computational framework to study problems of the form:

<p align="center">
  <img src="./doc/usage/img/gov_eq.png" alt="gov_eq.png" width=80% />
</p>

in a flexible and efficient manner, where <pre> <b>q</b> in R^n </pre> is a vector of ``n`` real-valued unknowns, ``t`` is time, and <pre> f(<b>q</b>) </pre> is a generic function of **``q``** which may include differential and algebraic operators.

The ability of dNami to clearly separate the problem statement from its numerical implementation (often a major time sink in research laboratories) is rooted in the flexibility of the Python language so as to let the user define her/his own system of balance laws in the most natural way (i.e. using a human-readable syntax), which is then interpreted in Fortran to build a computationally-efficient library of the equation above which is callable from Python. Users can then easily interact with their own system of balance laws, including at runtime, thereby making it possible to integrate solutions to the equation above with other tools and libraries (e.g. optimisation and stability tools) to fully explore the properties of the system, seamlessly from small to large-scale calculations.



## Dependencies

The following Python packages are needed:

* python3.x
* mpi4py
* numpy
* scons

In addition the following software is needed:

* MPI (OpenMPI or Intel-MPI recommended)
* A Fortran compiler gfortran or the Intel compiler ifort

Optionally, to run the provided plotting scripts to visualise the code output, matplotlib is required. 

## Quickstart guide

Check out the [documentation](https://dnami.readthedocs.io/en/latest/) to get a quickstart guide on running a case with dNami. 

## Test suite 

A Python script in the `exm/` folder called `test_all.py` will check the required dependencies (and will inform the user of any missing items) and then run through the list of available cases to test them. For each case, the script will copy the files to the appropriate location to generate the code, run the case and output values for a comparison with the reference values on the repository. Each case should output a `PASS` status once it is complete. If this is not the case, the `log.test` file generated during the testing process should clarify the reasons for the failure (e.g. missing dependence, lack of resources, ...). The testing script currently assumes that the user is testing the code on a machine with at least 4 available cores.  

## Generating the documentation locally  

The repository contains a documentation inside the **doc** directory.
In order to generate the html documentation the following Python packages are needed:
1. Sphinx
2. sphinx-rtd-theme
3. pydata-sphinx-theme
4. sphinxcontrib.bibtex

They can be installed using the following command:
```bash
pip3 install -U Sphinx
pip3 install -U sphinx-rtd-theme
pip3 install -U pydata-sphinx-theme
pip3 install -U sphinxcontrib.bibtex
```

Build the documentation by changing into the doc directory and executing the following command:
```bash
make html
```
After building the documentation the *_build/html* directory should contain the index.html startpage.


## Community guidelines 
-----------------------

### Contributing examples and test cases 

We look forward to seeing dNami used to solve problems from many different fields. If you would like to contribute an example, please prepare an `exm/*`-like folder with the equations, numerics and a compute.py as well as a short write-up of the problem and dNami results in the `Test cases and validation` section of the documentation. Then submit a pull request so that we can review it for acceptance. Thank you for helping us grow the example section.  

### Extending the documentation

If you want to add additional content to the documentation, add/edit *.rst* files in the *doc/usage*
directory and also update *doc/index.rst* (if necessary). 

### Issues and support  

If you discover a bug or an issue with the code, please open an issue with a clear write-up of the problem and steps to replicate it for us to investigate. 

## Authors
-------
See the AUTHORS file.

## License
-------
dNami is released under the New BSD License (see the LICENSE file for details).
