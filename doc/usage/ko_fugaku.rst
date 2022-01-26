The OIST ARM (Ko-Fugaku) cluster
********************************
This tutorial will explain how to run the 2D-Vortex-Advection example 
on the OIST ARM cluster. 
The OIST ARM cluster is currently not integrated into the SLURM job submission
system. For this reason jobs can only executed in interactive mode. 

* :ref:`system-overview` 


.. _system-overview:

System overview and login
=========================

The OIST ARM cluster consists of eight nodes, each with 48 cores.

   .. code-block:: bash

      kofugaku01 ~ kofugaku08

Login to Ko-Fugaku

   Any of the 8 nodes can be used to login. The following command
   can be used to login into the first node.

   .. code-block:: bash

      ssh -X YourName@kofugaku01.oist.jp

The following storage systems are available:

   +-------+----------+
   | Path  | Capacity |
   +=======+==========+
   | /work | 50TB     |
   +-------+----------+

   You can create your own working directory on **/work** and
   store/process your data directly on **/work**.


Running dNami in interactive mode
=================================

#. Clone dNami from the github repository

   In your home directory execute the following command

   .. code-block:: bash

    git clone git@github.com:oist/dNami.git
    git checkout oist_devNA_dNamiF_restored_nompi4py

   .. Warning:: Currently only the **oist_devNA_dNamiF_restored_nompi4py** branch is working on the ARM architecture.
      All other branches are missing the correct compiler flags inside the setup.scons.

#. Install the MPI4PY and scons Python packages

   .. code-block:: bash

    python3 -m pip install --user scons
    env MPICC=/opt/FJSVstclanga/v1.1.0/bin/mpifcc python -m pip install --user mpi4py

   .. Note:: There is no need to install Numpy because it is already installed.

#. Select the correct compiler inside the setup.scons file
   
   In order to compile dNami on the ARM architecture two compiler
   setups are available:

   +------------------------+------------------------------------------+
   | Flag inside setup.scons| Description                              |
   +========================+==========================================+
   | Ko-Fugaku-Fujitsu      | Uses the Fujitsu compiler to build dNami |
   +------------------------+------------------------------------------+
   | Ko-Fugaku-GCC          | Uses the GCC compiler to build dNami     |
   +------------------------+------------------------------------------+
  
   Change the beginning of the **setup.scons** file as shown below to select
   the Fujitsu compiler (remove **#** in front of the flag):

   .. code-block:: bash

    #COMPILER="GCC"
    #COMPILER="Intel"
    #COMPILER="Ko-Fugaku-GCC"
    COMPILER="Ko-Fugaku-Fujitsu"
    #COMPILER="Fujitsu"
    #COMPILER="FLANG"

   .. note:: 

    Using the Fujitsu compiler is recommended. 

#. Change into the **dNami/exm/2d_vortex_advection** directory, copy the two files genRhs.py and rhs.py to the src/generate directory

   .. code-block:: bash

      cp genRhs.py ../../src/generate
      cp rhs.py ../../src/generate

#. Change into the **src** directory and run the script

   .. code-block:: bash

      ./install_clean.sh

#. If your environment is setup correctly it should compile and build the dNami library. Add the dNami library to your path, from inside the src directory execute the command


   .. code-block:: bash

      source env_dNami.sh

#. Set the number of OpenMP threads to 1

   .. code-block:: bash

      export OMP_NUM_THREADS=1

#. Change to the **dNami/exm/2d_vortex_advection/** directory and run the example with the following command

   .. code-block:: bash

      mpirun -n 24 python3 compute.py

#. In case you want to run on multiple nodes you need to create a hostfile and pass it to the mpirun command:

   .. code-block:: bash

      mpirun -n 384 -hostfile hostfile -x PYTHONPATH python3 compute.py

   The hostfile will contain the names of the other nodes. You can copy the content of the hostfile shown below
   and save it inside a file named hostfile:

   .. code-block:: bash
      :caption: Content of the hostfile
      :name: hostfile

      kofugaku01 slots=48
      kofugaku02 slots=48
      kofugaku03 slots=48
      kofugaku04 slots=48
      kofugaku05 slots=48
      kofugaku06 slots=48
      kofugaku07 slots=48
      kofugaku08 slots=48

   In order to avoid interference of multiple users running jobs simultaneously on the cluster. Users can change 
   the hostfile to only include specific nodes. 
   For example, user1 uses the nodes from 1 ~ 4 and user2 from 5 ~ 8.
   In such a case the hostfile of user1 would look like this:

   .. code-block:: bash
      :caption: Content of the hostfile for user1
      :name: user1hostfile

      kofugaku01 slots=48
      kofugaku02 slots=48
      kofugaku03 slots=48
      kofugaku04 slots=48

   And the hostfile of user2 would look like this:

   .. code-block:: bash
      :caption: Content of the hostfile for user2
      :name: user2hostfile

      kofugaku05 slots=48
      kofugaku06 slots=48
      kofugaku07 slots=48
      kofugaku08 slots=48


