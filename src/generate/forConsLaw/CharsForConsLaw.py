import sympy as sym
import sys

sym.init_printing(use_latex=True,wrap_line=False)

def characteristics(eqn):

    #################################################################################
    #
    #    Some automative characteristic derivations for systemes of conservation laws
    #
    #################################################################################
    
    #####################################
    #  Deriving the quasi-linear form
    #####################################
    
    # Example of the Euler's equations in cartesian coordinates
    
    # General variables for the Euler equations:
    
    x,y,z,t     = sym.symbols(['x','y','z','t']   ,Real=True)
    xi,eta,zeta  = sym.symbols(['xi','eta','zeta'],Real=True)    

    rho,u  ,v  ,w  ,et ,E   = sym.symbols(['rho','u'  ,'v'  ,'w'  ,'et' ,'E'],Real=True)
    rhou ,rhov ,rhow ,rhoet = sym.symbols(['rhou' ,'rhov' ,'rhow' ,'rhoet'],Real=True)
    x_xi,y_xi,z_xi          = sym.symbols(['x_xi','y_xi','z_xi'],Real=True)    
    x_eta,y_eta,z_eta       = sym.symbols(['x_eta','y_eta','z_eta'],Real=True)  
    x_zeta,y_zeta,z_zeta    = sym.symbols(['x_zeta','y_zeta','z_zeta'],Real=True)  
    J = sym.Symbol('J',real=True)

    # Equation of states:
    #IG EOS
    
    gamma = sym.Symbol('gamma')
    p     = (gamma-1)*(rho*et-(rhou**2+rhov**2)/(2*rho))
    
    #Arbitrary EOS:
    # p   = sym.Function('p')(rho*et-0.5*(rhou**2+rhov**2)/rho,rho)
    # This general definition will introduce the Grüneisen index in the following derivations.
    
    # Defining the flux vector "Flx":
    #                  ---  
    #   ∂   ->    -->  ---     ->
    #   ──( f ) + Div( Flx ) = 0
    #   ∂t         
    # and the quasi-linear form and the "A matrix, A = [Ax,Ay,Az]""
    #
    #             -----------      
    #   ∂   ->    -----------  ->    ->
    #   ──( f ) + (A. \nabla)( f ) = 0
    #   ∂t         
    #
    

    # Flx_x for cartesian coordinates
    if eqn == 'Euler':
        Flx_x   = sym.Matrix([[rho*u],
                             [rho*u*u+p],
                             [rho*v*u],
                             [(rho*et+p)*u]])
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])
        Flx_y   = sym.Matrix([[rho*v],
                             [rho*u*v],
                             [rho*v*v+p],
                             [(rho*et+p)*v]])
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])    
        # Flx_x for curvilinear coordinates    
   
        Flx_x   = sym.Matrix([[1./J*rho*        ( u*x_xi + v*x_eta  )         ],
                              [1./J*(rho*u*     ( u*x_xi + v*x_eta  )+p*x_xi) ],
                              [1./J*(rho*v*     ( u*x_xi + v*x_eta  )+p*x_eta)],
                              [1./J*((rho*et+p)*( u*x_xi + v*x_eta  ))       ]])    
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])    
        
        Flx_y   = sym.Matrix([[1./J*rho*        ( u*y_xi+v*y_eta  )          ],
                              [1./J*(rho*u*     ( u*y_xi+v*y_eta  ) +p*y_xi )],
                              [1./J*(rho*v*     ( u*y_xi+v*y_eta  ) +p*y_eta)],
                              [1./J*((rho*et+p)*( u*y_xi+v*y_eta  ))        ]])    
                             # [rho*w],
                             # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]]) 
    
     #   Flx_z   = sym.Matrix([[rho* ( u*x_zeta + v*y_zeta  )],
     #                        [rho*u*( u*x_zeta + v*y_zeta  )+p*x_zeta],
     #                        [rho*v*( u*x_zeta + v*y_zeta  )+p*y_zeta],
     #                        [(rho*et+p)*( u*x_zeta + v*y_zeta  )]])    
     #                        # [rho*w],
     #                        # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])         

    else:
        print('Eqn not implemented yet.')
        import sys;sys.exit()                                          
        
    prim = [rho*u ,rho*v,u       ,v       ,rho*et]
    cons = [rhou  ,rhov ,rhou/rho,rhov/rho, rhoet] 
    switch = {}
    switchb= {}
    for p,c in zip(prim,cons):
        switch[p] = c
        switchb[c]= p
    
    
    Flx_x = Flx_x.subs(switch,simultaneous=True,evaluate=False)
    Flx_y = Flx_y.subs(switch,simultaneous=True,evaluate=False)    

    Q    = sym.Matrix([[rho],
                       [rhou],
                       [rhov],
                       # [w],
                       [rhoet]])
    
    # Ax, aY and Az are obtained from the Jacobian of the Flx_x, Flx_y and Flx_z compenant 
    Ax = Flx_x.jacobian(Q)
    Ay = Flx_y.jacobian(Q)     

    Ax = Ax.subs(switchb,simultaneous=True,evaluate=False) # to switch back to "velocities" expressions
    Ay = Ay.subs(switchb,simultaneous=True,evaluate=False) # to switch back to "velocities" expressions

    p = sym.Symbol('p')
    
    Q_CS = sym.Matrix([[rho],
                       [rho*u],
                       [rho*v],
                       [0.5*rho*(u*u+v*v)+1.0/(gamma-1)*p]])
                       # [rho*w],
                       # [0.5*rho*(u*u+v*v+w*w)+1.0/gamma_m1*P]])
    
    Q2    = sym.Matrix([[rho],
                       [u],
                       [v],
                       # [w],
                       [p]])
    
    
    M   = Q_CS.jacobian(Q2)
    Mm1 = M.inv() 

    Atild = Mm1*Ax*M
    Btild = Mm1*Ay*M    

    c = sym.Symbol('c')
    p     = (gamma-1)*(rho*et-(u**2+v**2)*rho/2)
    
    cexp = sym.simplify(gamma*p/rho)

    Atild = sym.simplify(Atild)
    Btild = sym.simplify(Btild)

    eq = c**2/gamma/p*rho
    
    etsub = sym.solve(eq-1,et)
    
    Atild = Atild.subs(et,etsub[0])
    Atild = sym.simplify(Atild)

    Btild = Btild.subs(et,etsub[0])
    Btild = sym.simplify(Btild)
    
    rho = sym.Function('rho')(x,y,z,t)
    u = sym.Function('u')(x,y,z,t)
    v = sym.Function('v')(x,y,z,t)
    w = sym.Function('w')(x,y,z,t)
    p = sym.Function('p')(x,y,z,t)
    
    Phi_x = sym.Matrix([rho,u,v,p]).diff(x)
    Phi_y = sym.Matrix([rho,u,v,p]).diff(y)
    
    Li,charSpeeds_x, Pleft_x   = charspeed_Lis(Atild,Phi_x)
    Mi,charSpeeds_y, Pleft_y = charspeed_Lis(Btild,Phi_y)

    # rho = sym.Function('rho')(xi,eta,zeta,t)
    # u = sym.Function('u')    (xi,eta,zeta,t)
    # v = sym.Function('v')    (xi,eta,zeta,t)
    # w = sym.Function('w')    (xi,eta,zeta,t)
    # p = sym.Function('p')    (xi,eta,zeta,t)    
    
    # Phi_xi = sym.Matrix([rho,u,v,p]).diff(xi)
    # Phi_eta = sym.Matrix([rho,u,v,p]).diff(eta)
    
    # Li,charSpeeds_xi, Pleft_xi   = charspeed_Lis(Atild,Phi_xi)
    # Mi,charSpeeds_eta, Pleft_eta = charspeed_Lis(Btild,Phi_eta)
 
    # sym.pprint(Li)
    
    return {'xi' :[Li,charSpeeds_x , Pleft_x],
            'eta':[Mi,charSpeeds_y, Pleft_y]}


#
#  Get characteristic speeds and eigenvectors
# 
def charspeed_Lis(Atild,Phi):
    (Sm1,D) = Atild.diagonalize()
    
    S = Sm1.inv()
    
    LeftEV = [ v[2] for v in Atild.left_eigenvects()]
    
    EV     = []
    for v in Atild.left_eigenvects():
      for i in range(0,v[1]):
        EV.append(v[0])
    
    
    Pleft = []
    for row in LeftEV:
      for a in row: 
        Pleft.append(a)             
    
    Pleft = sym.Matrix(Pleft)
    
    LiSym = Pleft*Phi
    
    LiFinal = sym.Matrix([e*l for e,l in zip(EV,LiSym)])

    return LiFinal,EV,Pleft

if __name__ == '__main__':
    # name,fbeg,fend,nstep = sys.argv
    characteristics('Euler')