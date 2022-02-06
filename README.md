# dNami
**dNami** [di:nɑ:mi:] is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes on structured meshes. 
dNami uses MPI and cache blocking techniques to speed up stencil-based operations.

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

## Community guidelines 
-----------------------

### Contributing 

We look forward to seeing dNami used to solve problems from many different domain. If you would like to contribute an example, please prepare an `exm/*`-like folder with the equations, numerics and a compute.py as well as a short write-up of the problem and dNami results in the `Test cases and validation` section of the documentation. Then submit a pull request so that we can review it for acceptance.   

### Issues and support  

If you discover a bug or an issue with the code, please open an issue with a clear write-up of the problem and steps to replicate it for us to investigate. 

## Authors
-------
See the AUTHORS file.

## License
-------
dNami is released under the New BSD License (see the LICENSE file for details).
