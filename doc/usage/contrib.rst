Contributing
***************************

We very much look forward to growing the dNami developer/user community. The dNami project is at its early stage and we have not demonstrated its full potential yet. However, it is mature enough to benefit from the momentum of a community of enthusiasts. We therefore welcome contributions. At the moment we are keen to grow the areas below.

Contributing examples and test cases
####################################

We look forward to seeing dNami used to solve problems from many different fields. If you would like to contribute an example, please prepare an ``exm/*``-like folder with the equations, numerics and a ``compute.py`` as well as a short write-up of the problem and dNami results in the :doc:`/usage/testcases/` section of the documentation. Then submit a pull request so that we can review it for acceptance. Thank you for helping us grow the example section.  

Extending the documentation
####################################

If you want to add additional content to the documentation, add/edit *.rst* files in the ``doc/usage``
directory and also update ``doc/index.rst`` (if necessary). Please make sure the documentation compiles locally before submitting a pull request. 

In order to generate the html documentation the following Python packages are needed:

	1. Sphinx
	2. sphinx-rtd-theme
	3. pydata-sphinx-theme
	4. sphinxcontrib.bibtex

They can be installed using the following command:

.. code-block:: bash

	pip3 install -U Sphinx
	pip3 install -U sphinx-rtd-theme
	pip3 install -U pydata-sphinx-theme
	pip3 install -U sphinxcontrib.bibtex


Build the documentation by changing into the ``doc`` directory and executing the following command:

.. code-block:: bash

	make html

After building the documentation the ``_build/html`` directory should contain the index.html startpage.


Issues and support
##################################

If you discover a bug or an issue with the code, please open an issue with a clear write-up of the problem and steps to replicate it for us to investigate. 
