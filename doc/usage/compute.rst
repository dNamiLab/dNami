Writing your own compute.py
***************************

The functions called in the `compute.py` are documented in the API reference section. However, this section aims to elucidate the order and reasoning behinds the various steps when setting up a case, allocating memory and passing information from the Python layer to the Fortran layer. 

Each of the following code blocks assume that the ``dnami`` Python library has been imported as ``dn``.

**Creating the tree** 

The first ``dnami`` library function called in the `compute.py` is: 

.. code-block:: python

        dtree = dn.create_tree()

This creates a dictionary of dictionaries each containing information about different aspects of the computation henceforth referred to as the 'tree'. Many dictionaries are initially empty but some contain information derived for the pseudo-code translation and compilation steps such as a list of the solved and stored variables or the number of halo points (determined by the largest finite-difference or filter stencil). For example, inspection of the `src/dnami.py` source shows an equation-related dictionary being populated with the list of solved variables:
 
.. code-block:: python

        for v in varsolved:
                dtree['eqns']['qvec']['solved'].append([v,varname[v]])   

where ``varsolved`` was imported from information generated from the `rhs.py` file. In the case of the 1D Euler equations solved in the quickstart guide, this list would contain `rho, u, et` which are the density, the velocity and the total energy.  After this tree creating step, the information about the computational parameters have to be give to the tree e.g. the grid size, the MPI domain decomposition, etc. For instance, the number of x-direction grid points ``nxgb`` can be set in the relevant section of the tree as:

.. code-block:: python

        dtree['grid']['size']['nxgb'] = nxgb

Note that this is the global number of points, not the processor specific number. The split across ``nxpr`` processors in the x-direction can be specified as:

.. code-block:: python

        dtree['mpi']['split']['nxpr'] = nxpr


**Initialising the Message Passing Interface**

If the total number of processors is not one, then the next step involves setting up the Message Passing Interface via: 

.. code-block:: python

        dtree = dn.start_mpi(dtree) 

With the information supplied in the previous section, the appropriate number of MPI processes are setup with the communicators for exchanging information between neighboring subdomains and specific ones for the boundary conditions if required. Currently, dNami uses MPI4PY. For more details the reader is referred to the `MPI4PY documentation <https://mpi4py.readthedocs.io/en/stable/>`_. The user can then obtain splitting-derived information such as the local number of points in the x-direction from the updated tree: 

.. code-block:: python

        nx = dtree['mpi']['dMpi'].nx 

Note that ``dtree['mpi']['dMpi']`` is a class.


**Allocating memory**

In this step, the memory used by each subprocess for the run parameters and the data (i.e. the solved and stored variables) is allocated:  

.. code-block:: python

        dtree = dn.allocate(dtree) 

Three main elements are allocated: a set of integer parameter (e.g. number of halo points, number of grid points, number of variables etc) which are used for memory reference purpose in the Fortran layer, a set of float parameters (e.g. grid spacing, time step, run constants, etc) and the data used and/or output during the run (see the difference between solved, stored and static variables). `Views <https://numpy.org/doc/stable/reference/generated/numpy.ndarray.view.html>`_ on these allocated memory regions are created so that the user can fill it (e.g. with the initial conditions) or perform operations with it (e.g. output the min/max of a given field). These views are added to the tree. The user can then create an alias to the views; for example, referring again to the 1D Euler case: 

.. code-block:: python 

        rho = dtree['eqns']['qvec']['views']['rho'] # density view
        u   = dtree['eqns']['qvec']['views']['u']   # velocity view
        et  = dtree['eqns']['qvec']['views']['et']  # total energy view

The user can then set the initial velocity field to zero:

.. code-block:: python 

       u[:] = np.float64(0.) 

To be clear, this operation does not create a new numpy array, it zeros the portion of the already-allocated memory that corresponds to the velocity variable. 


**Passing information to the Fortran layer**

A set of aliases for the three aforementioned arrays are created:

.. code-block:: python

        intparam,fltparam,data = (dtree['libs']['fort']['integers'],
                                  dtree['libs']['fort']['floats'],
                                  dtree['libs']['fort']['data'])

These memory references are then passed to the Fortran layer when calling the functions compiled with f2py e.g. when advancing the solution in time during the sub-RK steps:


.. code-block:: python

	dn.dnamiF.time_march(intparam,fltparam,data)  

The integer parameters (which are organised in a set pre-defined order) are used to read and modify the correct portion of the memory corresponding to ``data``. 


**Starting your own compute**

To create your own compute, we suggest that you start from an existing example that is closest to your desired case and tailor it to your needs. 

