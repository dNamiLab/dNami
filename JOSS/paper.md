---
title: '\texttt{dNami}: a framework for solving systems of balance laws using explicit numerical schemes on structured meshes.'
tags:
  - Python
  - Finite-difference 
  - Time-integration
  - Computational Physics
  - High Performance Computing
  - FORTRAN 
authors:
  - name: Nicolas Alferez^[co-first author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0003-3482-4529 
    affiliation: 1 
  - name: Emile Touber^[co-first author] ^[corresponding author] # note this makes a footnote saying 'co-first author'
    orcid: 0000-0002-7553-0351 
    affiliation: "2, 3"
  - name: Stephen D. Winn 
    orcid: 0000-0001-5972-4271 
    affiliation: "2, 3"
  - name: Yussuf Ali  
    affiliation: 2
affiliations:
 - name: Laboratoire DynFluid, Conservatoire National des Arts et Métiers, Paris, France  
   index: 1
 - name: Okinawa Institute of Science and Technology, Okinawa, Japan 
   index: 2
 - name: Department of Mechanical Engineering, Imperial College London, London, UK
   index: 3
date: 2 February 2022
bibliography: paper.bib

---

# Summary

Many physical systems obey or can be modelled by a system of balance laws taking the form of unsteady partial differential equations. Over the years, many research codes have been developed to solve such sets of equations, each requiring consequential development time despite many of them obeying a similar structure. \texttt{dNami} aims to provide a framework to rapidly set up and solve equations of the form 

\begin{equation} \label{eq:gov_eq}
	\frac{\partial\mathbf{q}}{\partial t} \textbf{q} = \textbf{RHS}(\mathbf{q})
\end{equation}

where $\textbf{q} \in \mathbb{R}^N$ with associated (time-dependant) boundary conditions on structured meshes while minimising problem-specific code development and still providing efficient and large-scale computing capabilities. 

\texttt{dNami} is aimed at both researchers who wish to rapidly set up a computation without having to worry too much about the coding aspects and those who wish to work with a high-performance code that they can expand and/or tailor to specific applications. For the former group, a syntax to symbolically define the governing equations and boundary conditions, an automatic symbolic-to-Fortran translation, a set of pre-implemented numerical methods all resulting in a native Python pre-compiled module generated automatically, provide a way of setting up computations and interacting conveniently with data at runtime without having to delve into Fortran code. For the latter group, \texttt{dNami} provides control over many performance aspects and easy integration of custom Fortran routines into the code generation step. Both groups can run efficient computations that can be started on a workstation and scaled up to a cluster.  

# State of the field

Code developments aiming at producing numerical solvers for problems given by \autoref{eq:gov_eq} typically follow two different strategies, most of the time involving two different communities. For simple enough types of \autoref{eq:gov_eq}, generic solutions provided by high-level languages may be used — notable examples are found in the vast offer provided by the scientific Python community. One of the reasons for the Python language success in computational science is its versatility to researchers' computational needs (e.g. diverse types of data structure, various interactions with data at runtime using the large offer from the community). A major drawback being its inability to tackle problems involving a large number of degrees of freedom, that typically requires highly efficient parallel capabilities. Note that solutions involving pre-compiled Python modules, see SciPy project [@SciPy2020], or just-in-time compilation, see Numba project [@lam2015numba] exist but are most of the time restricted to workstation workflow and incompatible with large scale computing of complex problems. The researcher targeting \autoref{eq:gov_eq} requiring large scale computing capabilities typically has to tackle the tedious problem of High-Performance Computing (HPC) development on his own. This is time consuming for a non-specialist and often results in conservative technical solutions that do not comply with the rapid evolution of hardware architecture that comes with continuous change of parallelisation paradigms. Another possibility is to collaborate with HPC specialists for the HPC-layer of the solver. Several examples of such fruitful joint efforts are available in the literature, see for instance @2021Streams in Computational Fluid Dynamics (CFD). However, such strategy often results in highly technical source codes, targetting at specific problem that are difficult to evolved in time and to understand for a non-HPC specialist. A relatively recent trends aiming at maintaining HPC capabilities while providing user and problem-specific flexibilities is using Domain-Specific-Languages (DSL) libraries developed by HPC specialists to tackle general problems like \autoref{eq:gov_eq}. Examples of such approaches can be found in CFD, a well-known highly HPC-resources demanding domain of computational physics [@2014PyFR; @2021HTR; @2021OpenSBLI]. Other DSL-based solvers target directly the \autoref{eq:gov_eq} [@2020Dedalus]. While DSL approach do provide versatility to physical problems and efficient adaptability to modern hardware architectures, they require users to learn the new DSL and drastically change paradigm in their developing approach. Besides, most DSL-based solutions rely on compiled binary executables produced from a low-level programming language (typically C or Fortran) to achieve performances. They do not provide the flexibility of pre-compiled Python modules from the user point of view (especially at runtime). \texttt{dNami} aims at conciliating the traditional DSL-based approach with all the advantages of a Python pre-complied module for targetting general problems like \autoref{eq:gov_eq}. 

# Features 

At the core of \texttt{dNami} is the translation of symbolic expressions written in high-level Python language to discretised equations in low-level Fortran language. \texttt{dNami} employs explicit schemes to discretise differential operators. For the temporal derivative of \autoref{eq:gov_eq}, a low-storage 3$^{rd}$ order Runge-Kutta scheme is used (other explicit schemes may easily be implemented). Spatial derivatives are discretised using finite differences of arbitrary orders provided by the user. A choice between standard and optimised schemes [such as those in @bogey2004family] is available. Both the governing equations in the form of \autoref{eq:gov_eq} and the boundary conditions can all be specified symbolically. \texttt{dNami} automatically deals with stencil-size and order reduction close to boundaries using the user-specified symbolic equations, which removes the need for time-consuming and often problem-specific code development. 

The source-to-source translation is performed by a set of Python functions using regular expressions (see `genKer.py`). The produced discretised version of \autoref{eq:gov_eq} is then inserted into appropriate do-loops included in Fortran template files by pre-processing techniques. This simple yet effective strategy makes it possible for the HPC-layer to be tailored at the template-file level independently of \autoref{eq:gov_eq}. Finally, the resulting Fortran source code is compiled as a shared library and optimised through the auto-optimisation process of modern Fortran compilers. A Python module is then created with a high-level interface of the shared library functions using F2PY [@2009f2py]. It is important to stress that researchers can still follow a traditional development workflow, in a pure Fortran environment, either at the shared library level or inside the high-level interface. 

\texttt{dNami} solves \autoref{eq:gov_eq} on structured grids. Parallelisation is then ensured via classical domain decomposition techniques through point-to-point MPI communications using  `mpi4py` [@2021mpi4py]. Efficient vectorisation of intense stencil-based computations are achieved via manual loop unrolling (yet automatically done by `genKer.py`) and cache-blocking techniques following recommendations from @2015andreolli. 

Python is used at runtime to set up the run parameters and initial conditions. More importantly the time loop is exposed to Python, giving the user the freedom to interact with the computation inside a Runge-Kutta step at run-time and to plug-in external libraries and/or output custom values with in-place data read and write between Python and Fortran. \texttt{dNami}'s Python interface thus allows easy integration of pre-processing, co-processing and post-processing tools. 

# Current \texttt{dNami} applications

\texttt{dNami} is currently being used to solve a wide variety of different physics problems such as ideal and non-ideal gasdynamics, (magneto-)hydrodynamic flows, shallow-water equations and global coupled air-water meteotsunami simulations. 


![Example computation of a global simulation of the water height variation due to the January 2022 Tonga volcano explosion.](earth_water.png){ width=80% }

# Acknowledgements

?

# References
