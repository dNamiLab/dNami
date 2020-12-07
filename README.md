![example workflow name](https://github.com/oist/dNami/workflows/1D_entropy_testcase/badge.svg)

# dNami

dNami is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes. 
dNami uses MPI, OpenMP and cache blocking to speed up the calculation.

## How-To generate the documentation
The repository contains a documentation inside the **doc** directory.
In order to generate the html documentation the following Python packages are needed:
1. Sphinx
2. sphinx-rtd-theme

They can be installed using the following command:
```bash
pip3 install -U Sphinx
pip3 install -U sphinx-rtd-theme
```

Build the documentation by changing into the doc directory and executing the following command:
```bash
make html
```
After building the documentation the *_build/html* directory should contain the index.html startpage.

### Extending the documentation

If you want to add additional content to the documentation, add/edit *.rst* files in the *doc/usage*
directory and also update *doc/index.rst* (if necessary). 

## Quickstart Guide

### Required software packages for running dNami
The following `Python` packages are needed:
1. `python3.x`
2. `mpi4py`
3. `matlibplot`
4. `numpy`
5. `scons` 

In addition the following software is needed:
1. MPI (`OpenMPI` or `Intel-MPI` recommended)
2. A Fortran compiler `gfortran` or the Intel compiler `ifort`

### Run the 2D-Vortex-Advection example
In order to verify your setup run the `2d_vortex_advection` example with
the following steps:

1. Change into the `dNami/exm/2d_vortex_advection` directory, copy the two files `genRhs.py` and `rhs.py` to the
`src/generate` directory:
     * `cp genRhs.py ../../src/generate/`
     * `cp rhs.py  ../../src/generate/`

2. Change into the `src` directory run the script: 
     * `./install_clean.sh`
     * If your environment is setup correctly it should compile and build the dNami library.

3. Add the dNami library to your path, from inside the `src` directory execute the command:
     * `source env_dNami.sh`

4. Change to the `dNami/exm/2d_vortex_advection/` directory and run the example with the
following command:
     * `mpirun -n 24 python3 compute.py` 


## Authors
-------
See the AUTHORS file.

## License
-------
dNami is released under the New BSD License (see the LICENSE file for details).
