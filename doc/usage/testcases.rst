Test cases and validation
*************************

One-dimensional cases
---------------------

1) Periodic entropy wave advection

This case solves the one-dimensional gas dynamics equations on a periodic domain:

.. math::

   \dfrac{\partial }{\partial t} \begin{pmatrix} \rho  \\ \rho u  \\ \rho e_t \end{pmatrix}  + \dfrac{\partial }{\partial x} \begin{pmatrix} \rho u   \\ \rho u^2 + p   \\ u ( \rho e_t + p) \end{pmatrix}   = 0

The files for this case are located in the ``exm/1d_entropywave`` folder. The initial non-dimensional density, pressure and velocity are set to :math:`(\rho/\rho_0, p/p_0, u/c_0) = (1,1,0.5)`. A sinusoidal perturbation of the density field is imposed. The case is run for 10 times the characterstic flow time. This test case compliments the quick-start guide in the sense that this case is periodic. No specification of a boundary condition in a given direction in the ``genRhs.py`` defaults to the domain being periodic in that direction. Periodicity is enforced via the ``swap`` operations. The resulting x-t diagram for the density field is shown in :numref:`xt_1d_per` for the 10 characteristic flow times.  



.. _xt_1d_per: 
.. figure:: img/xt_1d_periodic.png
   :width: 70%
   :align: center

   x-t diagram of the density perturbation for the one-dimensional Euler equations on a periodic domain. The blue lines indicate the flow speed. 

Two-dimensional cases
---------------------

1) Periodic vortex advection on a wavy mesh

This case solves the two-dimensional gasdynamics equations in curvilinear coordinates on a doubly-periodic domain using a wavy mesh for a weakly conservative formulation:

.. math::

   \dfrac{\partial }{\partial t} Q + \dfrac{1}{J} \left( \dfrac{\partial \eta}{\partial y} \dfrac{\partial }{\partial x} F_x - \dfrac{\partial \xi}{\partial y} \dfrac{\partial }{\partial x} F_y - \dfrac{\partial \eta}{\partial x} \dfrac{\partial }{\partial y} F_x + \dfrac{\partial \xi}{\partial x} \dfrac{\partial }{\partial y} F_y \right) = 0

where:

.. math::

   Q   = \begin{pmatrix} \rho  \\ \rho u \\ \rho v  \\ \rho e_t \end{pmatrix}, \ 
   F_x = \begin{pmatrix} \rho u \\ \rho u^2 + p \\ \rho u v \\ u ( \rho e_t + p) \end{pmatrix}, \
   F_y = \begin{pmatrix} \rho v  \\ \rho u v  \\ \rho v^2 +p    \\ v ( \rho e_t + p) \end{pmatrix}.

and :math:`J` is the Jacobian of the transformation between computational and physical space:

.. math::
   J \equiv \left( \dfrac{\partial \xi}{\partial x} \dfrac{\partial \eta}{\partial y} - \dfrac{\partial \eta}{\partial x}\dfrac{\partial \xi}{\partial y} \right)^{-1}

The computational space :math:`(x,y)` is related to physical space :math:`(\xi, \eta)` with the mapping: 

.. math::

   \xi  = \xi_0 + x L_x + A_x \sin( 2 \pi y) \\
   \eta = \eta_0 + y L_y + A_y \sin( 4 \pi x)

where the distances are specified relative to the vortex radius :math:`r_v`:

.. math::

   \xi_0 = -12 r_v, \ \eta_0 = -6 r_v, \ L_x = 24 r_v, \ L_y = 12 r_v, \  A_x   = 0.4 r_v, \ A_y = 1.6 r_v.

.. _2d_wmesh: 
.. figure:: img/2d_wavymesh_mesh.png
   :width: 70%
   :align: center

   Colormap showing the initial density distribution with vorticity contours shown in black. The mesh is displayed in light grey with :math:`(n_x,n_y)=(160,80)`. 


The metrics are computed with the same finite difference stencil and order as the derivatives in the governing equations. A vortex is initialised at the center of the domain at :math:`(\xi_c, \eta_c)=(0,0)`. The initial flow field is then specified as:

.. math::

   \left\{
   \begin{matrix}
   u(x,y,t=0) = u_0 \left( 1 - \dfrac{M_v}{M_i} \dfrac{\eta - \eta_c}{r_v} e^{(1-r^2)/2} \right) \\ 
   v(x,y,t=0) = v_0 \left( \dfrac{M_v}{M_i} \dfrac{\xi - \xi_c}{r_v} e^{(1-r^2)/2} \right) \\ 
   \end{matrix}
   \right.

The pressure and density are initialised based on isentropic ideal gas relations:

.. math::

   \left\{
   \begin{matrix}
   \rho(x,y,t=0) = \rho_0 \left( 1 - \dfrac{\gamma -1}{2} M_v^2 e^{(1-r^2)/2} \right)^{\dfrac{1}{\gamma - 1}} \\ 
   p(x,y,t=0)    = p_0    \left( 1 - \dfrac{\gamma -1}{2} M_v^2 e^{(1-r^2)/2} \right)^{\dfrac{\gamma}{\gamma - 1}} \\ 
   \end{matrix}
   \right.

The baseflow and vortex speed are specified via the Mach numbers :math:`M_i=0.5` and :math:`M_v=0.5` respectively. The mesh and initial density condition are shown in :numref:`2d_wmesh`. The case is run at a fixed grid size and timestep for various finite-difference schemes and orders. The same 11-point, 10 :sup:`th` order filter is used for every case. The results for standard finite difference schemes from 2 :sup:`nd` to 10 :sup:`th` order are shown in :numref:`2d_wmesh_ord`.  

.. _2d_wmesh_ord:
.. figure:: img/2d_wavymesh_FD.png
   :width: 50%
   :align: center

   Comparison of results for various finite-difference stencils and orders after :math:`t = 10 u_0/L_x` (i.e 10 vortex-travel times) for mesh size :math:`(n_x,n_y)=(160,80)`.      

2) Non-reflective vortex advection throught the boundaries

This case solves the two-dimensional advection of a vortex through the boundaries of the domain using a non-reflective characteristic-based boundary condition implementation. The two-dimensional Euler equations are:  

.. math::

   \dfrac{\partial }{\partial t} \begin{pmatrix} \rho  \\ \rho u \\ \rho v  \\ \rho e_t \end{pmatrix}  + \dfrac{\partial }{\partial x} \begin{pmatrix} \rho u   \\ \rho u^2 + p \\ \rho u v    \\ u ( \rho e_t + p) \end{pmatrix}  + \dfrac{\partial }{\partial y} \begin{pmatrix} \rho v   \\ \rho u v \\ \rho v^2 + p    \\ v ( \rho e_t + p) \end{pmatrix} = 0

supplemented with the ideal gas law:

.. math::
   
   p = \delta \rho \left[e_t - \dfrac{1}{2} ( u^2 + v^2) \right] 

The edge and corner boundaries are updated using a locally one-dimensional non-reflective boundary condition. For example, the upper boundaries are computed using the following expressions

.. math::

   \dfrac{\partial }{\partial t} 
   \left. \begin{pmatrix} \rho  \\ \rho u \\ \rho v  \\ \rho e_t \end{pmatrix} \right|_{x, y = L_y}  =  
    - \begin{pmatrix} d_1  \\ 
      u d_1 + \rho d_2 \\ 
      v d_1 + \rho d_3 \\ 
      (e_t + p/\rho + c_p/\alpha_v) d_1 + \rho u d_2 + \rho v d_3 + c_pd_4 /( \alpha_v  c^2) 
      \end{pmatrix}

where:

.. math::

    \begin{pmatrix} d_1  \\ d_2 \\ d_3  \\ d_4 \end{pmatrix} = 
    \begin{pmatrix} (\mathcal{L}_1 + \mathcal{L}_4 )/ c^2 + \mathcal{L}_2 \\ \mathcal{L}_3 \\ (\mathcal{L}_4 - \mathcal{L}_1) /(\rho c)  \\ \mathcal{L}_1 + \mathcal{L}_4 \end{pmatrix}, \ 
    \begin{pmatrix} \mathcal{L}_1  \\ \mathcal{L}_2 \\ \mathcal{L}_3   \\ \mathcal{L}_4 \end{pmatrix} =  
    \begin{pmatrix} \dfrac{1}{2} \max(v-c,0) \left( \dfrac{\partial p}{\partial y} - \rho c \dfrac{\partial v}{\partial y} \right)  \\ \max(v,0) \left( \dfrac{\partial \rho}{ \partial y} - \dfrac{1}{c^2} \dfrac{\partial p}{\partial y} \right) \\ \max(v,0) \dfrac{\partial u}{\partial x}   \\ \dfrac{1}{2} \max(v+c,0) \left( \dfrac{\partial p}{\partial y} + \rho c \dfrac{\partial v}{\partial y} \right)  \end{pmatrix}   


The initial flow field is set using:

.. math::

   \left\{
   \begin{matrix}
   \rho(x,y,t=0) = \rho_0, \\ 
   u   (x,y,t=0) = u_{0} - \dfrac{\partial \psi}{\partial y },  \\
   v   (x,y,t=0) = v_{0} + \dfrac{\partial \psi}{\partial x },  \\
   e_t (x,y,t=0) = (p_0 + p')/(\delta \rho_0) + \dfrac{1}{2} \left( u^2 + v^2 \right) 
   \end{matrix}
   \right.

where the derivatives of the potential :math:`\psi` and the pressure fluctuation :math:`p'` are set by:

.. math::

   \left\{
   \begin{matrix}
   & p'(x,y) = -\dfrac{\rho_0 \Gamma ^2}{2 R^2} e^{-r^2/R^2} , \\ 
   & \dfrac{\partial \psi}{\partial y }(x,y) = - \dfrac{y-y_0}{R^2} \Gamma e^{-r^2/(2R^2)},  \\
   & \dfrac{\partial \psi}{\partial x }(x,y) = - \dfrac{x-x_0}{R^2} \Gamma e^{-r^2/(2R^2)},  \\
   \end{matrix}
   \right.

The vortex is initially centered in the domain i.e. :math:`(x_0,y_0)=(0.5L_x, 0.5L_y)`


.. figure:: img/2d_vortexexit_drho.png
   :width: 100%

   Density fluctuations at various times during the interaction of the vortex with the non-reflective boundary. Vertical velocity contours are shown (with values in the range only at the start of the simulation).

Three-dimensional cases
-----------------------

1) Ideal gas Couette flow
