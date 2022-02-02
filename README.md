# dNami
**dNami** is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes. 
dNami uses MPI, OpenMP and cache blocking to speed up the calculation.

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

Check out the documentation [HERE] to get a quickstart guide on running a case with dNami. 

## How-To generate the documentation
The repository contains a documentation inside the **doc** directory.
In order to generate the html documentation the following Python packages are needed:
1. Sphinx
2. sphinx-rtd-theme
3. pydata-sphinx-theme

They can be installed using the following command:
```bash
pip3 install -U Sphinx
pip3 install -U sphinx-rtd-theme
pip3 install -U pydata-sphinx-theme
```

Build the documentation by changing into the doc directory and executing the following command:
```bash
make html
```
After building the documentation the *_build/html* directory should contain the index.html startpage.

### Extending the documentation

If you want to add additional content to the documentation, add/edit *.rst* files in the *doc/usage*
directory and also update *doc/index.rst* (if necessary). 

## Authors
-------
See the AUTHORS file.

## License
-------
dNami is released under the New BSD License (see the LICENSE file for details).
