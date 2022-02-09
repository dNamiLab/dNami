.. dNami documentation master file, created by
   sphinx-quickstart on Thu Oct 29 10:22:56 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dNami's documentation!
*********************************

**dNami** [di:nɑ:mi:] is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes on structured meshes. dNami uses MPI, OpenMP, loop-unrolling and cache-blocking techniques to speed up stencil-based operations.

.. only:: html

   .. _earth:
   .. figure:: usage/img/earth.gif 
      :width: 70%
      :align: center

      Animation of the atmospheric and water-height disturbance due to the Tonga volcano explosion in January 2022 computed with dNami. Check out a higher resolution version on YouTube: https://www.youtube.com/watch?v=sIivtcx1Al0

Motivation and scope
====================

The time evolution of a variety of physical and biological processes may be described by systems of balance laws, which, if given appropriate initial and boundary conditions, dictate the future states of the systems. For instance, systems of balance laws invoking mass, momentum and energy have been incredibly successful at providing meaningful insights to the future states of realistic systems in physics (e.g. fluid dynamics). Yet, experimenting numerically with such systems still requires much implementation time. dNami was created so that more research time is spent exploring the dynamical properties of the system of balance laws of interest to the user, and less time is wasted on its numerical implementation across the whole computational spectrum, from the initial small-scale exploratory work on a workstation to the final large-scale computations on national clusters. Thus, dNami is a computational framework to study problems of the form:

.. math::

   \begin{equation} \label{eq:gov_eq}
   \frac{\partial\textbf{q}}{\partial t} = \textbf{F}(\textbf{q}) \,\, + \,\, \mbox{initial/boundary conditions},
   \end{equation}

in a flexible and efficient manner, where :math:`\textbf{q} \in \mathbb{R}^n` is a vector of :math:`n` real-valued unknowns, :math:`t` is time, and :math:`\textbf{F}(\textbf{q})` is a generic function of :math:`\textbf{q}` which may include differential and algebraic operators.

The ability of Nami to clearly separate the problem statement from its numerical implementation (often a major time sink in research laboratories) is rooted in the flexibility of the Python language so as to let the user define her/his own system of balance laws in the most natural way (i.e. using a human-readable syntax), which is then interpreted in Fortran to build a computationally-efficient library of the  equation above which is callable from Python. Users can then easily interact with their own system of balance laws, including at runtime, thereby making it possible to integrate solutions to the equation above with other tools and libraries (e.g. optimisation and stability tools) to fully explore the properties of the system, seamlessly from small to large-scale calculations.

To get started with dNami, please check out the :doc:`Quickstart guide </usage/quickstart>` to get set up with dependencies and run your first dNami case. 


.. toctree::
   :caption: Content
   :maxdepth: 2

   usage/quickstart
   usage/testcases
   usage/genRhs   
   usage/compute
   usage/postpro
   usage/api
   usage/syntax
   usage/performance
   usage/references

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
