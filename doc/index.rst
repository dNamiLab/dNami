.. dNami documentation master file, created by
   sphinx-quickstart on Thu Oct 29 10:22:56 2020.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to dNami's documentation!
*********************************

**dNami** [di:na:mi] is an open-source multi-language (Python, Fortran, C) framework for solving systems of balance laws using explicit numerical schemes on structured meshes. dNami uses MPI, OpenMP, loop-unrolling and cache-blocking techniques to speed up stencil-based operations.

.. only:: html

   .. _earth:
   .. figure:: usage/img/earth.gif 
      :width: 70%
      :align: center

      Animation of the atmospheric and water-height disturbance due to the Tonga volcano explosion in January 2022 computed with dNami. 

To get started with dNami, please check out the Quickstart guide to get set up with dependencies and run your first dNami case. 


.. toctree::
   :caption: Content
   :maxdepth: 2

   usage/quickstart
   usage/testcases
   usage/genRhs   
   usage/compute
   usage/api
   usage/syntax
   usage/references

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
