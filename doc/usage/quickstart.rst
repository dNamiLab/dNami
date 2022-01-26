Quickstart Guide
****************

Required software packages for running dNami
The following Python packages are needed:

* python3.x
* mpi4py
* matlibplot
* numpy
* scons

In addition the following software is needed:

* MPI (OpenMPI or Intel-MPI recommended)
* A Fortran compiler gfortran or the Intel compiler ifort


Run the 2D-Vortex-Advection example
===================================
.. highlight:: sh

In order to verify your setup run the 2d_vortex_advection example with the following steps:

1. Change into the dNami/exm/2d_vortex_advection directory, copy the two files genRhs.py and rhs.py to the src/generate directory::

    cp genRhs.py ../../src/generate
    cp rhs.py ../../src/generate

2. Change into the src directory run the script::

    ./install_clean.sh

3. If your environment is setup correctly it should compile and build the dNami library. Add the dNami library to your path, from inside the src directory execute the command::

    source env_dNami.sh

4. Change to the dNami/exm/2d_vortex_advection/ directory and run the example with the following command::

    mpirun -n 24 python3 compute.py

