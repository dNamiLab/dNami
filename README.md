[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6720593.svg)](https://doi.org/10.5281/zenodo.6720593)
[![Documentation Status](https://readthedocs.org/projects/dnami/badge/?version=latest)](https://dnami.readthedocs.io/en/latest/?badge=latest)
# dNami
**dNami** [di:nɑ:mi:] is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes on structured meshes. 
dNami uses MPI, loop-unrolling and cache blocking techniques to speed up stencil-based operations. Spatial derivatives are constructed using a customisable finite-difference formulation. 

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

in a flexible and efficient manner, where **``q``** ``in R^n`` is a vector of ``n`` real-valued unknowns, ``t`` is time, and ``f(``**``q``**``)`` is a generic function of **``q``** which may include differential and algebraic operators.

The ability of dNami to clearly separate the problem statement from its numerical implementation (often a major time sink in research laboratories) is rooted in the flexibility of the Python language so as to let the user define her/his own system of balance laws in the most natural way (i.e. using a human-readable syntax), which is then interpreted in Fortran to build a computationally-efficient library of the equation above which is callable from Python. Users can then easily interact with their own system of balance laws, including at runtime, thereby making it possible to integrate solutions to the equation above with other tools and libraries (e.g. optimisation and stability tools) to fully explore the properties of the system, seamlessly from small to large-scale calculations.


## Why is dNami needed? 

Code developments aimed at producing numerical solvers for Equation (1) typically follow two different strategies which usually involve two different communities. 

For small-enough problems of the type of Equation (1), generic solutions provided by high-level languages may be used. Notable examples may be found in the vast offer provided by the scientific Python community. One of the reasons for the Python language success in computational science is its versatility to researchers’ computational needs (e.g. diverse object and data-type structures, various interactions with data at runtime using the large offer of tools from the community). However, a common major drawback is the inability to easily tackle problems involving a large number of degrees of freedom, which typically require highly efficient parallel capabilities. Note that solutions involving pre-compiled Python modules, see [SciPy project](https://doi.org/10.1038/s41592-019-0686-2), or just-in-time compilation, see [Numba project](https://doi.org/10.1145/2833157.2833162), exist but are most of the time restricted to workstation workflow and incompatible with large-scale computing of complex problems. 

Thus, researchers targeting Equation (1) and requiring large-scale computing capabilities typically have to tackle the tedious problem of High-Performance Computing (HPC) development on their own. This is time consuming for a non-specialist and often results in conservative technical solutions that do not comply with the rapid evolution of hardware architecture that comes with continuous change of parallelisation paradigms. Alternatively, such researchers can collaborate with HPC specialists for the HPC-layer of the solver. Several examples of such fruitful joint efforts are available in the literature, see for instance [STREAmS](https://doi.org/10.1016/j.cpc.2021.107906) and [cuIBM](https://doi.org/10.21105/joss.00301) in Computational Fluid Dynamics (CFD). However, this approach often results in highly technical source codes, targeting specific problems that are difficult to modify in time (e.g. new ``f(``**``q``**``)``, new choice of boundary conditions). 

A relatively recent trend aimed at maintaining HPC capabilities whilst providing user and problem-specific flexibilities is the use of Domain-Specific-Languages (DSL) libraries developed by HPC specialists to tackle  Equation (1). Examples of such approaches can be found in CFD, a notoriously HPC-intensive domain of computational physics e.g. [PyFR](https://doi.org/10.1016/j.cpc.2014.07.011), [HTR](https://doi.org/10.1016/j.cpc.2020.107733), [OpenSBLI](https://doi.org/10.1016/j.cpc.2021.108063). Other DSL-based solvers directly target Equation (1) e.g. [Dedalus](https://doi.org/10.1103/PhysRevResearch.2.023068) and [Coral](https://doi.org/10.21105/joss.02978). Although DSL approaches provide the versatility of the physical problem to solve and the efficient adaptability to modern hardware architectures, they do require users to learn the new DSL and drastically change paradigm in their developing approach. In addition, most DSL-based solutions rely on compiled binary executables produced from a low-level programming language (typically C or Fortran) to achieve their performance. They do not provide the flexibility of pre-compiled Python modules from the user point of view (especially at runtime). dNami aims at reconciling the DSL-based approach with all the advantages of a Python pre-complied module so as to solve general problems like Equation (1) on both small and large scales.

## What does dNami offer? 

At the core of dNami is the translation of symbolic expressions written in high-level Python language to discretise equations in low-level Fortran language. dNami employs explicit schemes to discretise differential operators. For the temporal derivative of Equation (1), a low-storage third order Runge--Kutta (RK) scheme is used (other explicit schemes may easily be implemented). Spatial derivatives are discretised using finite differences of arbitrary orders provided by the user. A choice between standard and optimised schemes (such as those in [Bogey and Bailly, 2004](https://doi.org/10.1016/j.jcp.2003.09.003)) is available. Both the governing equation in the form given above and the boundary conditions are specified symbolically. dNami automatically deals with stencil-size and order reduction close to boundaries using the user-specified symbolic equations, which removes the need for time-consuming and often problem-specific code development.

The source-to-source translation is performed by a set of Python functions using regular expressions (see `genKer.py`). The produced discretised version of Equation (1) is then inserted into appropriate do-loops included in Fortran template files by pre-processing techniques. This simple yet effective strategy makes it possible for the HPC-layer to be tailored at the template-file level independently of Equation (1). Finally, the resulting Fortran source code is compiled as a shared library and optimised through the auto-optimisation process of modern Fortran compilers. A Python module is then created with a high-level interface of the shared library functions using [F2PY](https://doi.org/10.1504/IJCSE.2009.029165). It is important to stress that researchers can still follow a traditional development workflow, in a pure Fortran environment, either at the shared library level or inside the high-level interface. 

dNami solves Equation (1) on structured grids. Parallelisation is then ensured via classical domain decomposition techniques through point-to-point MPI communications using  [mpi4py](https://doi/org/10.1109/MCSE.2021.3083216). Efficient vectorisation of intense stencil-based computations are achieved via manual loop unrolling (yet automatically done by `genKer.py`) and cache-blocking techniques following recommendations from [Andreolli, 2015](https://doi.org/10.1016/B978-0-12-802118-7.00023-6).

Python is used at runtime to set up the run parameters and initial conditions. More importantly the time loop is exposed to Python, giving the user the freedom to interact with the computation inside the RK steps at run-time and to plug-in external libraries and/or output custom values with in-place data read and write between Python and Fortran. dNami's Python interface thus allows easy integration of pre-processing, co-processing and post-processing tools. 

## Citing the code

If you use **dNami** for a scientific publication, please cite it using the [Zenodo reference](https://zenodo.org/record/6720593):

Bibtex:

```
@software{dNami,
	author       = {Alferez, Nicolas and Touber, Emile and Winn, Stephen and Ali, Yussuf},
	title        = {{dNami: a framework for solving systems of balance laws using explicit numerical schemes on structured meshes}},
	month        = jun,
	year         = 2022,
	publisher    = {Zenodo},
	version      = 2,
	doi          = {10.5281/zenodo.6720593},
	url          = {https://doi.org/10.5281/zenodo.6720593}
} 
```


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
After building the documentation the ``_build/html`` directory should contain the ``index.html`` startpage.


## Community guidelines 
-----------------------

### Contributing examples and test cases 

We look forward to seeing dNami used to solve problems from many different fields. If you would like to contribute an example, please prepare an `exm/*`-like folder with the equations, numerics and a compute.py as well as a short write-up of the problem and dNami results in the `Test cases and validation` section of the documentation. Then submit a pull request so that we can review it for acceptance. Thank you for helping us grow the example section.  

### Extending the documentation

If you want to add additional content to the documentation, add/edit ``.rst`` files in the ``doc/usage``
directory and also update ``doc/index.rst`` (if necessary). 

### Issues and support  

If you discover a bug or an issue with the code, please open an issue with a clear write-up of the problem and steps to replicate it for us to investigate. 

## Authors
-------
See the AUTHORS file.

## License
-------
dNami is released under the New BSD License (see the LICENSE file for details).
