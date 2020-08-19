import func_nz.values  as val
import func_nz.eos as eos
import numpy as np
from func_nz.integrate  import wrapper

def baseflow_model(x,y):
    """drhodx, dpdx and dMdx for a given x position. Input rho and M. Hospital's rule used at throat assumes the geometry is at least C1!"""

    # -- Unpack inputs
    rho,p,M  = y[0],y[1],y[2]
    one      = np.float64(1.)
    hlf      = np.float64(0.5)

    # -- Area functions defined in baseflow()
    if x == x_th : 
        # -- Throat condition        
        coeff  = np.sqrt( hlf/eos.Gamma(rho,p)*geom.d2A_fxdx2(x) /geom.A_fx(x) )
        drhodx = -rho*coeff
        dpdx   = -rho/val.Zc*eos.eos_sos(rho,p)**2*coeff
        dMdx   = -(one-eos.Gamma(rho,p)-one/M**2)*M*coeff
    else:
        # -- Everywhere else 
        coeff    = M**2/(M**2-one)*dA_fxdx(x)/A_fx(x)
        drhodx   = -rho*coeff
        dpdx     = -rho/val.Zc*eos.eos_sos(rho,p)**2*coeff
        dMdx     = -(one-eos.Gamma(rho,p)-one/M**2)*M*coeff

    return [drhodx,dpdx,dMdx]       

def baseflow(xin,Ain,dAdin,d2Ad2in,qin,xti =None):
    """Calculate the baseflow variables Q=(rho_bar,u_bar,p_bar) for a set of inlet conditions (rho1,p1,M1)"""

    # ----------------------------------------- Do the actual integration 

    # --- Set area functions
    global A_fx, dA_fxdx, d2A_fxdx2, x_th
    A_fx = Ain
    dA_fxdx = dAdin
    d2A_fxdx2 = d2Ad2in 
    x_th = xti

    # --- Convert to inlet values
    rhoi, pi, Mi = qin[0], qin[1], qin[2]

    #print("Computing perfect expansion for rho0, p0, M0 = [{}, {}, {}]".format(rhoi,pi,Mi))

    # --- Compute from entrance to exit  
    sol     = wrapper(baseflow_model,xin,qin)
    rho_bar = sol[:,0]
    p_bar   = sol[:,1]
    M_bar   = sol[:,2]
    
    # --- Nan check
    if np.isnan(M_bar).any() or np.isnan(rho_bar).any() or np.isnan(p_bar).any():
        print('NaN detected in the baseflow, exiting ...')
        exit()

    u_bar   = M_bar * eos.eos_sos(rho_bar,p_bar)

    # --- Package everything and output  
    nxx        = np.size(xin)
    Qbar       = np.zeros((nxx,3),dtype=np.float64)
    Qbar[:,0]  = rho_bar
    Qbar[:,1]  = u_bar
    Qbar[:,2]  = eos.eos_e(rho_bar,p_bar) + np.float64(0.5)*u_bar**2 

    return Qbar

def deriv_isentrope(x,y):
        #Compute the pressure values along the isentrope
        rho      = x
        p,M      = y[0],y[1]
        one              = np.float64(1.)
        c        = eos.eos_sos(rho,p)

        #Derivatives 
        dpdrho = c**2/val.Zc
        dMdrho = M/rho*(one-eos.Gamma(rho,p)-one/M**2)

        return [dpdrho,dMdrho]

def compute_thermo(rho1,p1,M1,dbg =0,manual=False):
    global rhot,pt

    one = np.float64(1.)
    nS  = 500 

    if M1 > 1.: 

        find_M = False
        factor = 1.5
        while not find_M:
            #Compute values along the isentrope until a sonic point is found:
            rho_end   = rho1*factor                           #this should be adjusted better
            rho_range = np.linspace(rho1,rho_end,nS)

            sol = wrapper(deriv_isentrope,rho_range,[p1,M1])
            ps      = sol[:,0]
            Ms      = sol[:,1]

            if np.amin(Ms) <= 1.  and np.amax(Ms) >= 1.:
                    find_M = True
            else:
                    factor *= 1.5

        from scipy.interpolate  import interp1d
        ft    = interp1d(Ms,rho_range)
        rhot  = ft(one)
        ft    = interp1d(rho_range,ps)
        pt    = ft(rhot)

    elif manual: 

        find_M = False
        factor = 1.5
        while not find_M:
            #Compute values along the isentrope until a sonic point is found:
            rho_end   = rho1*factor                           #this should be adjusted better
            rho_range = np.linspace(rho1,rho_end,nS)

            sol = wrapper(deriv_isentrope,rho_range,[p1,M1])
            ps      = sol[:,0]
            Ms      = sol[:,1]

            if np.amin(Ms) <= 1.  and np.amax(Ms) >= 1.:
                    find_M = True
            else:
                    factor *= 1.5

        from scipy.interpolate  import interp1d
        ft    = interp1d(Ms,rho_range)
        rhot  = ft(one)
        ft    = interp1d(rho_range,ps)
        pt    = ft(rhot)


    else:

        find_M = False
        factor = 2.
        while not find_M:
            #Compute values along the isentrope until a sonic point is found:
            rho_end   = rho1/factor                           #this should be adjusted better
            rho_range = np.linspace(rho1,rho_end,nS)

            sol = wrapper(deriv_isentrope,rho_range,[p1,M1])
            ps      = sol[:,0]
            Ms      = sol[:,1]

            if np.amin(Ms) <= 1.  and np.amax(Ms) >= 1.:
                    find_M = True
            else:
                    factor *= 2.

        from scipy.interpolate  import interp1d
        ft    = interp1d(Ms,rho_range)
        rhot  = ft(one)
        ft    = interp1d(rho_range,ps)
        pt    = ft(rhot)


    return rhot,pt

def find_AT(ri,pi,Mi):

    rhot, pt = compute_thermo(ri,pi,Mi)
    jstar = rhot * eos.eos_sos(rhot, pt)
    ci = eos.eos_sos(ri,pi)

    Ai_Astar = jstar / (ri * Mi * ci)

    return Ai_Astar 
