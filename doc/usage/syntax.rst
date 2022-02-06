dNami syntax
************

This section introduces the dNami syntax. By using the dNami syntax the 
user can easily translate differential equations into finite difference code.

First order derivative syntax
-----------------------------

The tables below show how to express derivatives in the dNami syntax in each of the three spatial directions. 

.. warning::

    When specifying derivatives, make sure to leave at least one white space after the ']_1x', ']_1y' or ']_1z' notation (this is so that the pseudo-code can be correctly understood and translated). For example:
        
         **WRONG**  ... [ f ]_1x+ 1.0_wp  ...

         **CORRECT**  ... [ f ]_1x + 1.0_wp  ...

.. table:: First derivative notation expressed in the dNami syntax

   +--------------------------------------+------------------+--------------------------------+
   |        Mathematical notation         |  dNami notation  |          Description           |
   +======================================+==================+================================+
   | .. math::                            |                  |                                |
   |                                      |                  |                                |
   |    \dfrac{\partial f}{\partial x}    |     [ f ]_1x     | First derivative in x direction|
   +--------------------------------------+------------------+--------------------------------+
   | .. math::                            |                  |                                |
   |                                      |                  |                                |
   |    \dfrac{\partial f}{\partial y}    |     [ f ]_1y     | First derivative in y direction|
   +--------------------------------------+------------------+--------------------------------+
   | .. math::                            |                  |                                |
   |                                      |                  |                                |
   |    \dfrac{\partial f}{\partial z}    |     [ f ]_1z     | First derivative in z direction|
   +--------------------------------------+------------------+--------------------------------+

Second order derivative syntax
------------------------------

To specify second order derivatives, two ways are currently possible. The user can directly specify a second derivative (discretised as a second derivative) or by taking the first derivative twice as detailed below. The two approaches are mathematically equivalent but will yields different results when discretised. The curly-bracket '}' symbol is used when taking a derivative inside another derivative. This approach can also be applied to cross-derivates.   

.. table:: Second derivative notation expressed in dNami syntax

   +---------------------------------------------------------------+------------------+----------------------------------------+
   |        Mathematical notation                                  |  dNami notation  |          Description                   |
   +===============================================================+==================+========================================+
   | .. math::                                                     |                  |                                        |
   |                                                               |                  |                                        |
   |    \dfrac{\partial^2 f}{\partial x^2}                         |     [ f ]_2xx    | Second derivative in x direction       |
   +---------------------------------------------------------------+------------------+----------------------------------------+
   | .. math::                                                     |                  |                                        |
   |                                                               |                  |                                        |
   |    \dfrac{\partial}{\partial x}\dfrac{\partial f}{\partial x} |                  | Double first derivative in x direction |
   |                                                               |    [ {f}_1x ]_1x |                                        |
   +---------------------------------------------------------------+------------------+----------------------------------------+
   | .. math::                                                     |                  |                                        |
   |                                                               |                  |                                        |
   |    \dfrac{\partial}{\partial y}\dfrac{\partial f}{\partial x} |     [ f ]_2xy    | Cross-derivative in x  and y directions|
   +---------------------------------------------------------------+------------------+----------------------------------------+
   | .. math::                                                     |                  |                                        |
   |                                                               |                  |                                        |
   |    \dfrac{\partial}{\partial y}\dfrac{\partial f}{\partial x} |                  | First derivative in x and y direction  |
   |                                                               |    [ {f}_1x ]_1y |                                        |
   +---------------------------------------------------------------+------------------+----------------------------------------+

Higher order derivative syntax
------------------------------

To generate higher-order derivatives, the current stategy involves storing an intermediate derivative and then taking the derivative of that stored variable. This is illustrated in the 1D KdV equations in :doc:`/usage/quickstart` where the third order derivative of the field :math:`u` is computed by computing and storing the second order derivative and then taking the first derivative of that stored field when specifying the RHS: 

.. code-block:: python

	varstored = { 'u_xx' : {'symb':' [u]_2xx ','ind':1, 'static': False } }
	...
	RHS = {'u' : ' epsilon * u * [ u ]_1x + mu * [ u_xx ]_1x ',}
