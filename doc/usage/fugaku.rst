The Fugaku cluster
******************
This tutorial will explain how to run the 2D-Vortex-Advection example 
on the Fugaku cluster.

* :ref:`system-overview` 


Installing the required Python packages
=======================================

#. Install the Python packages

   Fugaku uses two different node types Intel (login) and ARM (compute) nodes.
   The main challange is to compile the dNami code correctly on the login node.
   dNami requires numpy for compilation, if numpy is installed using pip3 on the
   login node the installtion will interfer with the installed numpy version on
   the compute node. This problem can be solved by installing a different version
   of Python on the login node.
   For this reason an environment for managing different Python versions is recommended.
   The following steps use the **asdf** software package to install different Python versions.

   Clone the current asdf version from the github repository:

   .. code-block:: bash

    git clone https://github.com/asdf-vm/asdf.git ~/.asdf --branch v0.8.1

   Add the following lines to your .bashrc:

   .. code-block:: bash

    . $HOME/.asdf/asdf.sh
    . $HOME/.asdf/completions/asdf.bash

   .. Warning::  After modifying the .bashrc logout and login again to Fugaku

   Install Python version 3.9.0 and scons+numpy using the following commands:

   .. code-block:: bash

    asdf plugin-add python
    asdf install python 3.9.0
    asdf global python 3.9.0
    python3.9 -m pip install scons --user
    python3.9 -m pip install numpy --user
    python3.9 -m pip install matplotlib --user

#. Clone dNami from the github repository

   In your home directory execute the following command

   .. code-block:: bash

    git clone git@github.com:oist/dNami.git
    git checkout oist_devNA_dNamiF_restored_nompi4py

   .. Warning:: Currently only the **oist_devNA_dNamiF_restored_nompi4py** branch is working on the ARM architecture.
      All other branches are missing the correct compiler flags inside the setup.scons.

#. Select the correct compiler inside the setup.scons file
   
   In order to compile dNami on the ARM architecture the Fujitsu compiler needs to be used
   Change the beginning of the **setup.scons** file as shown below to select
   the Fujitsu compiler (remove **#** in front of the flag):

   .. code-block:: bash

    #COMPILER="GCC"
    #COMPILER="Intel"
    #COMPILER="Ko-Fugaku-GCC"
    #COMPILER="Ko-Fugaku-Fujitsu"
    COMPILER="Fujitsu"
    #COMPILER="FLANG"

#. Change into the **dNami/exm/2d_vortex_advection** directory, copy the two files genRhs.py and rhs.py to the src/generate directory

   .. code-block:: bash

      cp genRhs.py ../../src/generate
      cp rhs.py ../../src/generate

#. Change into the **src** directory and run the script

   .. code-block:: bash

      ./install_clean.sh

#. If your environment is setup correctly it should compile and build the dNami library. 

   Copy the following code into a file named **submit.sh** and place it into the **dNami/exm/2d_vortex_advection** directory.

   .. code-block:: bash

     #!/bin/sh -x
     #PJM -L  "node=1"                          # Number of assign node 8 (1 dimention format)
     #PJM -L  "rscgrp=small"                    # Specify resource group
     #PJM -L  "elapse=00:15:00"                 # Elapsed time limit 1 hour
     #PJM --mpi "max-proc-per-node=48"          # Maximum number of MPI processes created per node
     #PJM --mpi "proc=48"
     #PJM -s                                    # Statistical information output
     
     module purge
     module load Python3-CN/3.6.8
     module load lang/tcsds-1.2.33
     
     export OMP_NUM_THREADS=1
     
     export CURRENT=$PWD
     cd ../../src
     source ./env_dNami.sh
     cd $CURRENT
     
     #export PLE_MPI_STD_EMPTYFILE=off
     export LD_PRELOAD=/usr/lib/FJSVtcs/ple/lib64/libpmix.so
     #export FLIB_CNTL_BARRIER_ERR=FALSE
     
     rm *.err
     rm *.out
     rm *.log
     
     mpirun -stdout output.log -stderr error.log -n 48 python3 ./compute.py

   .. Warning:: A single Fugaku node has 48 cores, adjust the domain decomposition inside the computer.py.

#. Submit the job script

   Use the following command to submit the job submission script:

   .. code-block:: bash

     pjsub submit.sh

   Use the follwoing command to check the status of the job submission:

   .. code-block:: bash

     pjstat
