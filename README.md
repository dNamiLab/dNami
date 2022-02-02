# dNami
**dNami** is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes. 
dNami uses MPI and cache blocking to speed up the calculation.

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

## Test suite 

A Python script in the `exm/` folder called `test_all.py` will run through the list of available cases to test them. For each case, the script will copy the files to the appropriate location to generate the code, run the case and output values for a comparison with the reference values on the repository. Each case should output a `PASS` status once it is complete. If this is not the case, the `log.test` file generated during the testing process should clarify the reasons for the failure (e.g. missing dependence, lack of resources, ...). The testing script currently assumes that the user is testing the code on a machine with at least 4 available cores.  

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
