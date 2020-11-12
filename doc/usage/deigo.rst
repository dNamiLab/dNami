Howto run dNami on the OIST Deigo cluster in interactive mode
*************************************************************
This tutorial will explain how to run the 2D-Vortex-Advection example 
on the OIST cluster Deigo. 

.. highlight:: sh


1. Login to Deigo::

    ssh -X YourName@deigo.oist.jp

2. Generate a SSH key (Press the Enter key so select the default values)::

    ssh-keygen -t ed25519 -C "your_email@oist.jp"

3. Register the SSH key in our Github page	
    .. image:: img/g1.jpg
       :width: 100%

4. Active your new ssh key for OIST	
    .. image:: img/g2.jpg
       :width: 100%

5. Clone dNami from the github repository, in your home directory execute the following command::

    git@github.com:oist/dNami.git

6. Load Python version 3.7 with the following command::
    
    module load python/3.7.3

7. Install the Python make system scons::

    python3 -m pip install --user scons

8. Change into interactive mode (if the cluster is very busy it may take some time until your request will be executed)::

    srun -t 0-1 -p compute -C zen2 -c 20  --mem=16G --pty 

9. From the command in step 8 you should get some output similar to the output below::

    srun: job 3783215 queued and waiting for resources
    srun: job 3783215 has been allocated resources

10. Change into the dNami/exm/2d_vortex_advection directory, copy the two files genRhs.py and rhs.py to the src/generate directory::

     cp genRhs.py ../../src/generate
     cp rhs.py ../../src/generate

11. Change into the src directory run the script::

     ./install_clean.sh

12. If your environment is setup correctly it should compile and build the dNami library. Add the dNami library to your path, from inside the src directory execute the command::

     source env_dNami.sh

13. Set the number of OpenMP threads to 1::

     export OMP_NUM_THREADS=1

14. Change to the dNami/exm/2d_vortex_advection/ directory and run the example with the following command::

     mpirun --oversubscribe -n 24 python3 compute.py

15. The output can be visualized by using the live_view.py script. You can login to Deigo with a second terminal window (keep the first terminal open to run the code). For running dNami you **must** be in **interactive mode**, for visualizing the output you don't need to be in interactive mode (running live_view from a login node is ok).
You can distinguish between the two modes by looking at your bash::

     your_name@deigo-login1 src]$    "login" indicates that you are on a login node
     your_name@deigo011706 src]$     "deigo011706" indicates that you are in the interactive mode (instead of 011706 it could also be a different number)

16. Copy the live_view.py file to the example directory (assuming you are inside the directory dNami/exm/2d_vortex_advection) ::

     mkdir out/liv
     cp ../../pst/liv/live_view.py ./out/liv

17. Run live_view.py with the following command (from inside the out/liv directory)::

     python3 live_view.py

18. If no new window opens on your MacOS screen, you may need to install XQuartz: https://www.xquartz.org/index.html


