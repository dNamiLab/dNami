---
title: '\texttt{dNami}: a framework for solving systems of balance laws using explicit numerical schemes on structured meshes'
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

The time evolution of a variety of physical and biological processes may be described by systems of balance laws, which, if given appropriate initial and boundary conditions, dictate the future states of the systems. For instance, systems of balance laws invoking mass, momentum and energy have been incredibly successful at providing meaningful insights to the future states of realistic systems in physics (e.g.\ fluid dynamics). Yet, experimenting numerically with such systems still requires much implementation time. \texttt{dNami} (di:na:mi:) was created so that more research time is spent exploring the dynamical properties of the system of balance laws of interest to the user, and less time is wasted on its numerical implementation across the whole computational spectrum, from the initial small-scale exploratory work on a workstation to the final large-scale computations on national clusters. Thus, \texttt{dNami} is a computational framework to study problems of the form:

\begin{equation} \label{eq:gov_eq}
\frac{\partial\textbf{q}}{\partial t} = \textbf{RHS}(\textbf{q}) \,\, + \,\, \mbox{initial/boundary conditions},
\end{equation}

in a flexible and efficient manner, where $\textbf{q} \in \mathbb{R}^n$ is a vector of $n$ real-valued unknowns, $t$ is time, and $\textbf{RHS}(\textbf{q})$ is a generic function of $\textbf{q}$ which may include differential and algebraic operators.

The ability of \texttt{dNami} to clearly separate the problem statement from its numerical implementation (often a major time sink in research laboratories) is rooted in the flexibility of the Python language so as to let the user define her/his own system of balance laws in the most natural way (i.e.\ using a human-readable syntax), which is then interpreted in Fortran to build a computationally-efficient library of \autoref{eq:gov_eq} which is callable from Python. Users can then easily interact with their own system of balance laws, including at runtime, thereby making it possible to integrate solutions to \autoref{eq:gov_eq} with other tools and libraries (e.g.\ optimisation and stability tools) to fully explore the properties of the system, seamlessly from small to large-scale calculations.

# State of the field

Code developments aimed at producing numerical solvers for \autoref{eq:gov_eq} typically follow two different strategies which usually involve two different communities. 

For small-enough problems of the type of \autoref{eq:gov_eq}, generic solutions provided by high-level languages may be used. Notable examples may be found in the vast offer provided by the scientific Python community. One of the reasons for the Python language success in computational science is its versatility to researchers’ computational needs (e.g.\ diverse object and data-type structures, various interactions with data at runtime using the large offer of tools from the community). However, a common major drawback is the inability to easily tackle problems involving a large number of degrees of freedom, which typically require highly efficient parallel capabilities. Note that solutions involving pre-compiled Python modules, see SciPy project [@SciPy2020], or just-in-time compilation, see Numba project [@lam2015numba] exist but are most of the time restricted to workstation workflow and incompatible with large-scale computing of complex problems. 

Thus, researchers targeting \autoref{eq:gov_eq} and requiring large-scale computing capabilities typically have to tackle the tedious problem of High-Performance Computing (HPC) development on their own. This is time consuming for a non-specialist and often results in conservative technical solutions that do not comply with the rapid evolution of hardware architecture that comes with continuous change of parallelisation paradigms. Alternatively, such researchers can collaborate with HPC specialists for the HPC-layer of the solver. Several examples of such fruitful joint efforts are available in the literature, see for instance @2021Streams and @Krishnan2017 in Computational Fluid Dynamics (CFD). However, this approach often results in highly technical source codes, targeting specific problems that are difficult to modify in time (e.g.\ new $\textbf{RHS}(\textbf{q})$, new choice of boundary conditions). 

A relatively recent trend aimed at maintaining HPC capabilities whilst providing user and problem-specific flexibilities is the use of Domain-Specific-Languages (DSL) libraries developed by HPC specialists to tackle  \autoref{eq:gov_eq}. Examples of such approaches can be found in CFD, a notoriously HPC-intensive domain of computational physics [@2014PyFR; @2021HTR; @2021OpenSBLI]. Other DSL-based solvers directly target \autoref{eq:gov_eq} [@2020Dedalus; @Miquel2021]. Although DSL approaches provide the versatility of the physical problem to solve and the efficient adaptability to modern hardware architectures, they do require users to learn the new DSL and drastically change paradigm in their developing approach. In addition, most DSL-based solutions rely on compiled binary executables produced from a low-level programming language (typically C or Fortran) to achieve their performance. They do not provide the flexibility of pre-compiled Python modules from the user point of view (especially at runtime). \texttt{dNami} aims at conciliating the DSL-based approach with all the advantages of a Python pre-complied module so as to solve general problems like \autoref{eq:gov_eq} on both small and large scales.

# Features 

At the core of \texttt{dNami} is the translation of symbolic expressions written in high-level Python language to discretise equations in low-level Fortran language. \texttt{dNami} employs explicit schemes to discretise differential operators. For the temporal derivative of \autoref{eq:gov_eq}, a low-storage 3$^{rd}$ order Runge--Kutta (RK) scheme is used (other explicit schemes may easily be implemented). Spatial derivatives are discretised using finite differences of arbitrary orders provided by the user. A choice between standard and optimised schemes [such as those in @bogey2004family] is available. Both the governing equation in the form of "$\partial\textbf{q}/\partial t = \textbf{RHS}(\textbf{q})$" and the boundary conditions are specified symbolically. \texttt{dNami} automatically deals with stencil-size and order reduction close to boundaries using the user-specified symbolic equations, which removes the need for time-consuming and often problem-specific code development.

The source-to-source translation is performed by a set of Python functions using regular expressions (see `genKer.py`). The produced discretised version of \autoref{eq:gov_eq} is then inserted into appropriate do-loops included in Fortran template files by pre-processing techniques. This simple yet effective strategy makes it possible for the HPC-layer to be tailored at the template-file level independently of \autoref{eq:gov_eq}. Finally, the resulting Fortran source code is compiled as a shared library and optimised through the auto-optimisation process of modern Fortran compilers. A Python module is then created with a high-level interface of the shared library functions using F2PY [@2009f2py]. It is important to stress that researchers can still follow a traditional development workflow, in a pure Fortran environment, either at the shared library level or inside the high-level interface. 

\texttt{dNami} solves \autoref{eq:gov_eq} on structured grids. Parallelisation is then ensured via classical domain decomposition techniques through point-to-point MPI communications using  `mpi4py` [@2021mpi4py]. Efficient vectorisation of intense stencil-based computations are achieved via manual loop unrolling (yet automatically done by `genKer.py`) and cache-blocking techniques following recommendations from @2015andreolli.

Python is used at runtime to set up the run parameters and initial conditions. More importantly the time loop is exposed to Python, giving the user the freedom to interact with the computation inside the RK steps at run-time and to plug-in external libraries and/or output custom values with in-place data read and write between Python and Fortran. \texttt{dNami}'s Python interface thus allows easy integration of pre-processing, co-processing and post-processing tools. 

# Current \texttt{dNami} applications

\texttt{dNami} is currently being used by researchers to study linear and nonlinear wave phenomena in various systems, from engineering to biology and geo/astro-physical flows. For example, we solve the Navier--Stokes equations to study shock waves and compressible turbulence, reaction-diffusion equations to study dissipative solitons in biological networks, the magneto-hydrodynamic equations to study shock-entropy interactions in plasmas, and the shallow-water equations to study geophysical flows (e.g.\ internal gravity waves, tsunamis and meteotsunamis, two-dimensional compressible turbulence).

To illustrate the ability of \texttt{dNami} to solve new systems of balance laws on supercomputers in record times, we used the eruption of the Tonga volcano on 12 January 2022, which made tsunami warnings fail in Japan [@tonga]. A new system of balance laws to couple gravity-driven waves (shallow-water equations) in the ocean and the primary Lamb wave (compressible Euler equations) in the atmosphere following the eruption was created and computed on a global scale using available bathymetry data. The arrival times of the meteotsunami waves are found to be in good agreement with tide gages around the globe. A snapshot of the resulting animation is shown in \autoref{fig:tonga}. Starting from scratch the whole process was completed in a day owing to the flexibility with which \texttt{dNami} could handle the new set of equations, the choice of spherical coordinates, and importing geotiff files for the bathymetry data.

![Example computation of a global simulation of the atmopsheric and water-height disturbance due to the January 2022 Tonga volcano explosion. Time is shown in UT starting on January 15th. \label{fig:tonga}](tonga.png){ width=100% }

# Acknowledgements

We are grateful for the help and support provided by the Scientific Computing and Data Analysis section of Research Support Division at OIST. This work used computational resources of the supercomputer Fugaku provided by RIKEN through the HPCI System Research Project (Project ID: hp200198 and hp210186).

# References
