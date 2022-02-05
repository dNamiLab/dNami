Writing your own genRhs.py
***************************

The functions called in the ``genRhs.py`` are documented in the API reference section. However, this section aims to elucidate the order and reasoning behinds the various steps when setting up the Partial Differential Equations (PDEs) system to be solved together with Boundary Conditions (BCs) and specifying numerical parameters associated to the automatic discretisation of derivative operators. A section aiming at describing how loop-distribution and -fusion can conveniently be done without digging into the Fortran layer is also provided.

How to generate the discretised equations
########################################################

``genRhs.py`` is the user-level module that specifies the system of equations to be marched in time. For convenience, this process is often divided into two steps: 1) the PDEs and BCs specifications through a series of dictionaries and lists, and 2) the Fortran source code generation through appropriate calls to ``genKer.py`` user-level functions.

Assume that ones want to march in time the 2D compressible Euler equations given by:

.. math::

   \dfrac{\partial }{\partial t} \begin{pmatrix} \rho  \\ \rho u \\ \rho v  \\ \rho e_t \end{pmatrix}  + \dfrac{\partial }{\partial x} \begin{pmatrix} \rho u   \\ \rho u^2 + p \\ \rho u v    \\ u ( \rho e_t + p) \end{pmatrix}  + \dfrac{\partial }{\partial y} \begin{pmatrix} \rho v   \\ \rho u v \\ \rho v^2 + p    \\ v ( \rho e_t + p) \end{pmatrix} = 0

**Setting equations in rhs.py**
 
PDEs and BCs specifications are typically done in a separate file, hereafter named ``rhs.py``. This file should contain the set of PDEs provided through dictionaries (see below) and a series of specific lists needed by dNami. The first parameter to be put in ``rhs.py`` is the number of spatial dimensions:

.. code-block:: python

			dim = 2 # can be 1,2 or 3

*Dictionaries of equations*

dNami does not contain any preset symbols or limited forms of PDEs. It is the user responsibility to declare the list of symbols needed to write the PDEs that will be marched in time and the corresponding equations. This is done through the ``varname`` dictionary (needed by the framework):

.. code-block:: python

	varname  = {'rho' : 1,
		    'u'   : 2,
		    'v'   : 3,
		    'et'  : 4, 
		    }

The order set with the indexes in ``varname`` corresponds to data location in memory. Variables solved in the left-hand side are specified in ``varsolved``.

.. code-block:: python

	varsolved = ['rho','u','v','et']

Setting RHS equations is done through a dictionary of the form 

.. code-block:: python

	{'Variable Name': 'Pseudo-Code Eqns'}

The ``Variable Name`` must match one of the ``varname`` keys. ``Pseudo-Code Eqns`` refers to symbolic expressions using ``varname`` keys or Fortran compatible syntax (Fortran internal functions are allowed here). 

.. warning::

	Although dNami marches in times systems of the form:

	.. math::

   		\dfrac{\partial \textbf{q} }{\partial t} = \textbf{RHS}\left( \textbf{q} \right)

	The ``'Pseudo-Code Eqns'`` are written on the left in::

		{'Variable Name': 'Pseudo-Code Eqns'}

For the x derivative of the RHS in the 2D Euler equations we would write:

.. code-block:: python
	
	divFx = {'rho' : ' [ rho*u           ]_1x ', 
    		 'u'   : ' [ rho*u*u + p     ]_1x ', 
    		 'v'   : ' [ rho*v*u         ]_1x ', 
    		 'et'  : ' [ (rho*et + p )*u ]_1x ', 
    }

In this expression the pressure is introduced through a new symbol, ``'p'``, not defined in ``varname``. Two possibilities are offered by dNami in such cases. The first one is to provide an equation that relates ``'p'`` with ``varname`` variables, this is done through the ``varloc`` dictionary:

.. code-block:: python

        varloc = { 'e' : ' (et - 0.5_wp*u*u) ',                        
                   'p' : '       rho*e       ',                        
                   }

dNami will automatically replace any occurrence of ``'p'`` with the corresponding combination of ``varname`` variables in all treatment of ``'Pseudo-Code Eqns'`` provided to the kernel (through ``append_Rhs`` or ``genBC``).
Another option is to allocate static memory for ``'p'`` and compute ``'p'`` before filling the RHS, where only memory access to that location are done. This is done through the ``varstored`` dictionary:

.. code-block:: python
	
	varloc = { 'e' : ' (et - 0.5_wp*u*u) '}                      
	varstored = {'p' : {'symb': 'rho*e', 'ind':1 , 'static': True}

In this example, an equation is provided to compute ``'e'`` from ``varname`` and ``'p'`` is stored at the first location of the stored-data memory.



*List of solver parameters*

**The** ``append_Rhs`` **function**

[




**Compulsory steps**

.. code-block:: python

			from genKer import rhsinfo, genrk3, genrk3update, genFilter, genBC, append_Rhs, genbcsrc
			import os 
			
			wp = 'float64'


**Optional steps**

*Adding explicit filtering*

*Adding boundary conditions*

Advanced use: control of the Fortran loop distribution
######################################################


[TO DO Needs to comment on the following successive steps]

.. code-block:: python

	varloc       = {'e': '(et-0.5_wp*(u*u+v*v))'}
	varstored    = {'p'  : {'symb' : 'gamma_m1*rho*(e)', 'ind' : 1, 'static' : False}}
	SRC = {
   			'rho' :  ' [ rho*u          ]_1x ',
   			'u'   :  ' [ rho*u*u  + p   ]_1x ',
   			'v'   :  ' [ rho*u*v        ]_1x ',
   			'et'  :  ' [ (rho*et + p)*u ]_1x '
   			}

with 

.. code-block:: python

	append_Rhs(SRC, 3, 2, rhsname,vnamesrc_divFconv,update=False,rhs=rhs)

producing the following Fortran code (we could detail vnamesrc_divFconv and rhsname)

.. code-block:: fortran


	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	! [rho*u*u+p]_1x
	!
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	d1_FluXconv_dx_0_im1jk = q(i-1,j,indvars(1))*q(i-1,j,indvars(2))*q(i-1,j,indvars(2))+qst(i-1,j,indvarsst(1))
	
	d1_FluXconv_dx_0_ip1jk = q(i+1,j,indvars(1))*q(i+1,j,indvars(2))*q(i+1,j,indvars(2))+qst(i+1,j,indvarsst(1))
	
	d1_FluXconv_dx_0_ijk = -&
	          0.5_wp*d1_FluXconv_dx_0_im1jk+&
	          0.5_wp*d1_FluXconv_dx_0_ip1jk
	
	d1_FluXconv_dx_0_ijk = d1_FluXconv_dx_0_ijk*param_float(1)
	
	
	
	!***********************************************************
	!                                                           
	! Update RHS terms for d(rho u)/dt *************************
	!                                                           
	!***********************************************************
	
	
	rhs(i,j,indvars(2)) =   -  ( d1_FluXconv_dx_0_ijk ) 

and 

.. code-block:: fortran

	!***********************************************************
	!                                                           
	! Start computing stored variables *************************
	!                                                           
	!***********************************************************
	
	
	     do j=idloop(3),idloop(4) 
	 
	      do i=idloop(1),idloop(2) 
	
	
	!***********************************************************
	!                                                           
	! building stored variable p *******************************
	!                                                           
	!***********************************************************
	
	
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	!
	! gamma_m1*rho*(e)
	!
	!~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	
	
	
	!***********************************************************
	!                                                           
	! Update stored variables gamma_m1*rho*(e) *****************
	!                                                           
	!***********************************************************
	
	
	qst(i,j,indvarsst(1)) =  param_float(3 + 5)*q(i,j,indvars(1))*((q(i,j,indvars(4))-&
	                    0.5_wp*(q(i,j,indvars(2))*q(i,j,indvars(2))+&
	                    q(i,j,indvars(3))*q(i,j,indvars(3)))))
	
	     enddo
	   enddo