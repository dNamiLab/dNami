Tools
*************************

A certain number of tools are provided alongside the core of dNami. They are briefly presented in this section.

I/O Tools
---------------------

In the **pst/utils** folder, a **post_io.py** script is provided which contains three I/O functions:

* *load_ax()* which takes the path of the axes as an argument and returns the number of grid points and the axes arrays   
* *read_restart()* which takes the path of a core restart files and return the simulation time and the *q* array 
* *read_restart_wshells()* which takes the path of a core restart files and return the simulation time and the *q* array with the shell information if they exist

The location of this script is added to the python path when the dNami environment is sourced (in **src/env_dNami.sh**). The aforementioned functions can then be imported as:

.. code-block:: python

        from post_io import read_restart,read_restart_wshell, load_ax

