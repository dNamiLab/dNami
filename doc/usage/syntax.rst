The dNami syntax
****************

This section introduces the dNami syntax. By using the dNami syntax the 
user can easily translate differential equations into finite difference code.

Expressing differential equations in dNami
------------------------------------------

The tables below show how to express derivatives in the dNami syntax.

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

.. table:: Second derivative notation expressed in dNami syntax

   +----------------------------------------+------------------+---------------------------------+
   |        Mathematical notation           |  dNami notation  |          Description            |
   +========================================+==================+=================================+
   | .. math::                              |                  |                                 |
   |                                        |                  |                                 |
   |    \dfrac{\partial^2 f}{\partial x^2}  |     [ f ]_2xx    | Second derivative in x direction|
   +----------------------------------------+------------------+---------------------------------+
   | .. math::                              |                  |                                 |
   |                                        |                  |                                 |
   |    \dfrac{\partial^2 f}{\partial y^2}  |     [ f ]_2yy    | Second derivative in y direction|
   +----------------------------------------+------------------+---------------------------------+
   | .. math::                              |                  |                                 |
   |                                        |                  |                                 |
   |    \dfrac{\partial^2 f}{\partial z^2}  |     [ f ]_2zz    | Second derivative in z direction|
   +----------------------------------------+------------------+---------------------------------+

