---
title: '\texttt{dNami}: a framework for solving systems of balance laws using explicit numerical schemes.'
tags:
  - Python
  - finite-difference 
  - time-integration
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
 - name: Laboratoire DynFluid, Arts et Métiers ParisTech, Paris, France  
   index: 1
 - name: Okinawa Institute of Science and Technology, Okinawa, Japan 
   index: 2
 - name: Department of Mechanical Engineering, Imperial College London, London, UK
   index: 3
date: 2 February 2022
bibliography: paper.bib

---

# Summary

Many physical systems obey or can be modelled by a system of balance laws. Over the years, many research codes have been developed to solve such sets of equations, each requiring consequential development time despite many of them obeying a similar structure. \texttt{dNami} aims to provides a framework to rapidly set up and solve equations of the form 

\begin{equation} \label{eq:gov_eq}
	\frac{\partial}{\partial t} \textbf{q} = \textbf{RHS}(\textbf{q})
\end{equation}

where $\textbf{q} \in \mathbb{R}^N$ with associated (time-dependant) boundary conditions while minimising problem-specific code development.

\texttt{dNami} is aimed at both researchers who wish to rapidly set up a computation without having to worry too much about the coding aspects and those who wish to work with a high-performance code that they can expand and/or tailor to specific applications. For the former group, a syntax to symbolically define the governing equations and boundary conditions, an automatic symbolic-to-Fortran translation, a set of pre-implemented numerical methods and a Python interface provide a way of setting up computations without having to delve into Fortran code. For the latter group, \texttt{dNami} provides control over many performance aspects and easy integration of custom Fortran routines into the code generation step. Both groups can run efficient computations that can be started on a workstation and scaled up to a cluster.  


# Features 

At the core of \texttt{dNami} is the translation of symbolic expressions to discretised equations. Both the governing equations in the form of \autoref{eq:gov_eq} and the boundary conditions can all be specified symbolically. \texttt{dNami} automatically deals with stencil/order reduction close to boundaries using the user-specified symbolic equations which removes the need for time consuming and often problem-specific code development. 

For numerical discretisation, \texttt{dNami} employs an explicit finite-difference approach where a choice between standard and optimised schemes can be made (such as those in @bogey2004family). Users can also easily supply their own schemes and filter coefficients. Time integration is performed with an explicit low-storage 3$^{rd}$ order Runge-Kutta scheme. The symbolic equations are automatically translated to Fortran code with the user-specified numerical parameters.  

For performance, \texttt{dNami} uses MPI to efficiently make use of available computational resources. Cache-blocking [TODO]

A Python interface is used to set up the run parameters and initial conditions. This interface, which wraps around call to the Fortran layer, also allows the user to interact with the computation at run-time to plug-in external libraries and/or output custom values. \texttt{dNami}'s Python interface allows easy integration of pre-processing, co-processing and post-processing tools. 

# Current \texttt{dNami} applications

\texttt{dNami} is currently being used to solve a wide variety of different physics problems such as ideal and non-ideal gasdynamics, (magneto-)hydrodynamic flows, shallow-water equations and global coupled air-water meteotsunami simulations. 


![Example computation of a global simulation of the water height variation due to the January 2022 Tonga volcano explosion.](earth_water.png){ width=80% }

# Acknowledgements

?

# References
