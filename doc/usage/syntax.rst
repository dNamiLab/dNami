dNami syntax
************

This section introduces the dNami syntax. By using the dNami syntax the 
user can easily translate differential equations into finite difference code.

Expressing differential equations in dNami
------------------------------------------

The tables below show how to express derivatives in the dNami syntax in each of the three spatial directions. 

.. warning::

    When specifying derivatives, make sure to leave at least one white space after the ']_1x', ']_1y' or ']_1z' notation (this is so that the pseudo-code can be correctly understood and translated).

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

To specify second order derivatives, two ways are currently possible. The user can directly specify a second derivative (discretised as a second derivative) or by taking the first derivative twice as detailed below. The two approaches are mathematically equivalent but will yields different results when discretised. The curly-bracket '}' symbol is used when taking a derivative inside another derivative.   

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

To generate higher-order derivatives, the current stategy involves storing an intermediate derivative and then taking the derivative of that stored variable. 
