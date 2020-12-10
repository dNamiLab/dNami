How-To run dNami on the OIST Deigo cluster in interactive mode
**************************************************************
This tutorial will explain how to run the 2D-Vortex-Advection example 
on the OIST cluster Deigo. 

* :ref:`req-label` 
* :ref:`interactive-label`
* :ref:`batch-label`

.. _req-label:

Pre-requisites
==============

The following steps assume that you already have a Github account.
If you don't have one please go to www.github.com and create an account.

.. highlight:: sh


#. Login to Deigo

   .. code-block:: bash

      ssh -X YourName@deigo.oist.jp

#. Generate a SSH key 

   The next command may prompt some questions where to save the generated key, simply
   press the **Enter** key to select the default values.

   .. code-block:: bash

      ssh-keygen -t ed25519 -C "your_email@oist.jp"

   The next command will print your generated public key. The whole output of the
   command is needed in the next step when you add your key to your Github account.

   .. code-block:: bash

      cat ~/.ssh/id_ed25519.pub

#. Register the SSH key in our Github page	

   .. image:: img/g1.jpg
      :width: 100%

#. Activate your new ssh key for OIST	

   .. image:: img/g2.jpg
      :width: 100%

.. _interactive-label:

Running dNami from in interactive mode
======================================

#. Clone dNami from the github repository

   In your home directory execute the following command

   .. code-block:: bash

    git clone git@github.com:oist/dNami.git

#. Load Python version 3.7 with the following command
    
   .. code-block:: bash

    module load python/3.7.3

#. Install the Python make system scons

   .. code-block:: bash

    python3 -m pip install --user scons

#. Change into interactive mode
   
   If the cluster is very busy it may take some time until your request will be executed.
   Try the following command first.

   .. code-block:: bash
      :caption: 1
      :name: Try-1
  
      srun -t 0-1 -p short --ntasks 20  --mem=16G --pty bash -l
    
   If your request was successful you should see that your terminal prompt changed as shown below.
   Instead of **deigo-login*** it will show something similar to **deigo011706**

   .. code-block:: bash

      your_name@deigo-login1 ~]$  "login" indicates that you are on a login node
      your_name@deigo011706  ~]$  "deigo011706" indicates that you are in the interactive mode (instead of 011706 it could also be a different number)

   .. Caution:: Double check that you are on a compute node (interactive mode), running heavy workloads on login nodes is forbidden and
      may have an impact on the usage of the Deigo system.

   It may happen that the command in :ref:`Try-1` takes some time to be evaluated.
   You may also see some output similar to the output below. 

   .. code-block:: bash
      :caption: The change into the interactive mode was successful

      srun: job 3783215 queued and waiting for resources
      srun: job 3783215 has been allocated resources

   If the command in :ref:`Try-1` takes a long time, cancel the request by pressing Ctrl+c and try
   the following:

   .. code-block:: bash
      :caption: Try to change into the interactive mode on another partition
      :name: Try-2

      srun -t 0-1 -p compute --ntasks 20  --mem=16G --pty bash -l

   Some additional background information on the options and the available partitions on Deigo.
   (You can skip this for the moment)

   +------------+------------------------------------------------------+
   | Option     | Explanation                                          |
   +============+======================================================+
   | -t 0-1     | You want to use Deigo for **0** days and **1** hour  |
   +------------+------------------------------------------------------+
   | -p short   | You want to use the **short** partition              |
   +------------+------------------------------------------------------+
   | -C zen2    | You want to use the AMD CPUS                         |
   +------------+------------------------------------------------------+
   | --ntasks 20| You want to use 20 MPI processes                     |
   +------------+------------------------------------------------------+
   | --mem=16G  | Reserve 16 GB of RAM                                 |
   +------------+------------------------------------------------------+

   The following image shows the Deigo cluster partition layout, as a student you hava access
   to the **short** and **compute** partition.
   
   .. image:: img/deigo_overview.png
      :width: 45%
   .. image:: img/deigo_partition.png
      :width: 50%

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

#. The output can be visualized by using the live_view.py script. 

   You can login to Deigo with a second terminal window (keep the first terminal open to run the code). 
   For running dNami you **must** be in **interactive mode**, for visualizing the output you don't need to be in interactive mode (running live_view from a login node is ok).
   You can distinguish between the two modes by looking at your terminal prompt:

   .. code-block:: bash

      your_name@deigo-login1 ~]$   "login" indicates that you are on a login node
      your_name@deigo011706  ~]$   "deigo011706" indicates that you are in the interactive mode (instead of 011706 it could also be a different number)

#. Copy the live_view.py file to the example directory (assuming you are inside the directory **dNami/exm/2d_vortex_advection**)

   .. code-block:: bash

      mkdir out/liv
      cp ../../pst/liv/live_view.py ./out/liv

#. Run live_view.py with the following command (from inside the out/liv directory)

   .. code-block:: bash

     python3 live_view.py
   
   If no new window opens on your MacOS screen, you may need to install XQuartz: https://www.xquartz.org/index.html


#. You can exit the interactive mode by the following command

   .. code-block:: bash

      exit

   After exiting the interactive mode you are back on the Deigo login node

   .. code-block:: bash

      your_name@deigo011706  ~]$ exit
      your_name@deigo-login1 ~]$  



.. _batch-label:

Running dNami from a batch script
=================================

dNami can also be executed using a job submission script.
The script :ref:`batch-1`  can be used as a template, it assumes that **gFortran** and
**OpenMPI** are used to compile the program. 
Copy the code from below and save it in the same directory as your compute.py, use the
filename **deigo_script.sh**

   .. code-block:: bash
      :caption: Job script template
      :name: batch-1

      #!/bin/bash
      #SBATCH --job-name=YourJobName
      #SBATCH --partition=short
      #SBATCH -C zen2
      #SBATCH --time=01:20:00
      #SBATCH --mem=500G
      #SBATCH --ntasks=1024
      #SBATCH --cpus-per-task=1
      
      #SBATCH --threads-per-core=1
      #SBATCH --sockets-per-node=2
      #SBATCH --cores-per-socket=64
      #SBATCH --ntasks-per-node=128
      #SBATCH --ntasks-per-socket=64
      #SBATCH --ntasks-per-core=1
      #SBATCH --exclusive
      
      module purge
      module load python/3.7.3
      
      cd ../../src/
      source env_dNami.sh
      cd $SLURM_SUBMIT_DIR
      
      export OMP_NUM_THREADS=1
      
      srun --mpi=pmix python3.7 compute.py > output.log 2>&1

You can adjust the following options to match your settings in the compute.py. 

   +-----------+--------------------------------------------------------+
   | Option    | Explanation                                            |
   +===========+========================================================+
   | --job-name| Set a job name, if you run multiple jobs you can       |
   |           | distinguish between different jobs.                    |
   +-----------+--------------------------------------------------------+
   | --time    | You can set a rough (or precise) estimate how long your|
   |           | job may take to run. The more precise you are the      |
   |           | higher the propability that your job starts earlier.   |
   +-----------+--------------------------------------------------------+
   | --mem     | Set the amount of memory you want to use, this setting |
   |           | is per node. 500G is the maximum.                      |
   +-----------+--------------------------------------------------------+
   | --ntasks  | Set the number of MPI processes, this setting must     |
   |           | match the product of your *with_proc* setting inside   |
   |           | your compute.py.                                       |
   +-----------+--------------------------------------------------------+

Submit the job script from the same directory where you placed yout compute.py.

   .. code-block:: bash

     sbatch deigo_script.sh

You can check the status of your jobs by using the following command:

   .. code-block:: bash

     squeue

The job will write all the output to the file **output.log**. 