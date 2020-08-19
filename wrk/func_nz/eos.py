#All the variations on the equation of state 
import func_nz.values as val
import numpy  as np

# ----------------------------------------------------------------------------------------

gas_mdl   = val.gas_type

#--------- Ideal gas

if gas_mdl == 0 :

        def init_eos():
                val.Zc = np.float64(1.)
                return 

        def eos_p(rho,e):
                p = (val.gamma-1.)*rho*e
                return p 

        def eos_e(rho,p):
                e = p/(rho*(val.gamma-1.))
                return e

        def eos_sos(rho,p):
                c =  np.sqrt(val.gamma*p/rho)
                return c

        def eos_rho(p,T):
                rho = p/T
                return rho

        def eos_T(rho,p):
                T = p/rho
                return T

        def eos_p_T(rho,T):
                #Done
                p = rho*T
                return p 

        def eos_sos_T(rho,T):
                #Done 
                p = eos_p_T(rho,T)
                c = np.sqrt(val.gamma*p/rho)
                return c                


#--------- Van der Waals gas

elif gas_mdl == 1 :

        def init_eos():
                val.Zc = np.float64(3.)/np.float64(8.)
                return 

        def eos_T(rho,p):
                thr = np.float64(3.)
                egh = np.float64(8.)
                T   = (p+thr*rho**2)/egh/rho*(thr-rho)
                return T

        def eos_T_E(rho,e):
                thr    = np.float64(3.)
                egh    = np.float64(8.)
                n_o_e  = np.float64(9.)/np.float64(8.)
                T = val.delta*(e+n_o_e*rho)
                return T                

        def eos_p(rho,e):
                thr = np.float64(3.)
                egh = np.float64(8.)
                T   = eos_T_E(rho,e)
                p   = egh*rho*T/(thr-rho) - thr*rho**2
                return p 

        def eos_rho(p,T):
                roots = np.roots([3.0, -9.0, (8.0*T + p), -3.0*p])
                r1  = min(roots.real)
                for i in range(0,3):
                        if roots[i].imag  == 0.0 and roots[i].real > r1 :
                                r1 = roots[i].real
                return r1

        def eos_p_T(rho,T):
                thr = np.float64(3.)
                egh = np.float64(8.)
                p   = egh*rho*T/(thr-rho) - thr*rho**2
                return p 

        def eos_e(rho,p):
                n_o_e = np.float64(9.)/np.float64(8.)
                T         = eos_T(rho,p)
                e         = T / val.delta - n_o_e*rho
                return e

        def eos_sos(rho,p):
                one       = np.float64(1.)
                two       = np.float64(2.)
                thr       = np.float64(3.)
                n_o_e     = np.float64(9.)/np.float64(8.)
                T                 = eos_T(rho,p)
                c_squared = T*(one + val.delta)*( one / ( one - rho / thr) )**2 - two*rho*n_o_e
                c         = np.sqrt( c_squared )
                #if np.isnan(c).any():
                #        print('c^2 is negative ... in eos_sos  ')
                #        exit()
                return c

        def eos_sos_T(rho,T):
                one       = np.float64(1.)
                two       = np.float64(2.)
                thr       = np.float64(3.)
                n_o_e     = np.float64(9.)/np.float64(8.)
                c_squared = T*(one+val.delta)*( one / ( one - rho /thr) )**2 - two*rho*n_o_e
                c         = np.sqrt( c_squared )
                if np.isnan(c).any():
                        print('c^2 is negative ... in eos_sos_T ')
                        exit()          
                return c                

else: 
        print('Please select an equation of state (values.py -> gas_mdl == 0,1,2)')


# ---------------------------------------------------------------------------------------- 

# ------- Functions based on EoS

# ------------------------------------------------------------------------------------#
#Calculating Hugoniot (shock adiabat) lines
def hugoniot(rho1,p1,rho2):

        from scipy.optimize import fsolve

        #Internal energy at point 1 :
        e1 = eos_e(rho1,p1) 
        v1 = 1./rho1

        #Create p2 vector to be filled 
        p2   = np.zeros(rho2.size)
        idx  = 0
        pini = p1 

        for rhoh in rho2:

                vh = 1./rhoh

                #Create the function to zero
                def H_line(pg):

                        #Internal energy at point 2
                        e2 = eos_e(rhoh,pg)

                        #Hugoniot function to zero 
                        f  = e2 - e1 - 0.5*(pg+p1)*(v1-vh)*val.Zc

                        return f

                #Find corresponding p2 and save it      
                p2[idx] = fsolve(H_line,pini,xtol=1e-10)
                pini    = p2[idx]               
                idx     = idx + 1

        # print('----------------------------------')

        return p2                       

# ------------------------------------------------------------------------------------#
#Calculating Hugoniot (shock adiabat) lines
def shock_jump(rho1,p1,M1,dbg=0,exp=False):

        from scipy.optimize import fsolve

        e1   = eos_e  (rho1,p1) 
        c1   = eos_sos(rho1,p1)
        v1       = 1./rho1

        if dbg==1:
                print('Debug shock_jump [1]     ')
                print(e1,c1,rho1)

        def equations(w):
                rho2, p2 = w

                e2 = eos_e(rho2,p2) 
                v2 = 1./rho2
                j  = M1*c1/v1

                return (e2 - e1 - 0.5*(p2+p1)*(v1-v2)*val.Zc, p2-p1 + 1.0/val.Zc*j**2*(v2-v1) )

        #NOTE: this assumes a compression shock ... 
        if exp:
            wini=(rho1/2.5, p1/2.5)
        else:
            wini=(10.*rho1, 10.*p1)
                
        rho2, p2 =  fsolve(equations,wini,xtol=1e-10)   

        if dbg==1:
                print('Debug shock_jump [2]')
                print(rho2, p2)

        M2 = rho1*M1*c1/rho2/eos_sos(rho2,p2)
        if M2 > 1.:

            print('r2 / p2 / M2= {} {} {}'.format(rho2,p2,M2))
            print('Post shock Mach number greater than 1 ...')
            exit()

        return (rho2,p2)

# -------- Thermodynamic derivatives: 

def dcdrho(rho,p,drho = 1e-6,ig_debug=0):
        "Compute dcdrho at fixed p "
        c_right = eos_sos(rho+drho,p)
        c_left  = eos_sos(rho-drho,p)
        dcdrho  = (c_right - c_left)/(np.float64(2.)*drho)
        #Exact expression for an ideal gas              
        if ig_debug == 1:
                dcdrho = -0.5*eos_sos(rho,p)/rho
        return dcdrho

def dcdp(rho,p,dp = 1e-6,ig_debug=0):
        "Compute dcdp at fixed rho "
        c_right = eos_sos(rho,p+dp)
        c_left  = eos_sos(rho,p-dp)
        dcdp    = (c_right - c_left)/(np.float64(2.)*dp)
        #Exact expression for an ideal gas      
        if ig_debug == 1:
                dcdp = 0.5*eos_sos(rho,p)/p
        return dcdp

def dp2dr2f(p1,r1,r2,dbg = 1):
        """Returns the derivatives of p2 with respect to v2 along the Hugoniot line.
        Done using central difference. Adjust the deltav as necessary."""
        from scipy.interpolate  import interp1d

        #Compute the Hugoniot line up to and a bit further from v2
        delr  = 1e-6
        rvec  = np.linspace(r1,r2+np.float64(2.)*delr,val.nH)
        pvec  = hugoniot(r1,p1,rvec)
        pf    = interp1d(rvec,pvec)

        #Central difference 
        dp2dr2 = (pf(r2+delr)-pf(r2-delr))/(np.float64(2.)*delr)

        # if dbg == 1:
        #       print('ph(rho2) in dp2dr2=',pf(r2))

        return dp2dr2

def dp2dr1f(p1,r1,r2):
        """Returns the derivatives of p2 with respect to r1 along the Hugoniot line.
        Done using central difference. Adjust the delta as necessary."""
        from scipy.interpolate  import interp1d

        delr   = 1e-6

        rvec  = np.linspace(r1-delr,r2,val.nH)
        pvec  = hugoniot(r1-delr,p1,rvec)
        pfm   = interp1d(rvec,pvec)

        rvec  = np.linspace(r1+delr,r2,val.nH)
        pvec  = hugoniot(r1+delr,p1,rvec)
        pfp   = interp1d(rvec,pvec)

        #Central difference 
        dp2dr1 = (pfp(r2) - pfm(r2))/(np.float64(2.)*delr)

        return dp2dr1   

def dp2dp1f(p1,r1,r2):
        """Returns the derivatives of p2 with respect to p1 along the Hugoniot line.
        Done using central difference. Adjust the delta as necessary."""
        from scipy.interpolate  import interp1d

        delp   = 1e-6

        rvec  = np.linspace(r1,r2,val.nH)
        pvec  = hugoniot(r1,p1-delp,rvec)
        pfm   = interp1d(rvec,pvec)

        rvec  = np.linspace(r1,r2,val.nH)
        pvec  = hugoniot(r1,p1+delp,rvec)
        pfp   = interp1d(rvec,pvec)

        #Central difference 
        dp2dp1 = (pfp(r2) - pfm(r2))/(np.float64(2.)*delp)

        return dp2dp1   

def dedrho(rho,p):
        "Compute dedrho at fixed p "
        drho    = 1e-6
        e_right = eos_e(rho+drho,p)
        e_left  = eos_e(rho-drho,p)
        dedrho  = (e_right - e_left)/(np.float64(2.)*drho)
        return dedrho

def dedp(rho,p,dp = 1e-6,ig_debug=0):
        "Compute dedp at fixed rho "
        e_right = eos_e(rho,p+dp)
        e_left  = eos_e(rho,p-dp)
        dedp    = (e_right - e_left)/(np.float64(2.)*dp)
        return dedp

def dc_Tdrho(rho,T):
        drho    = 1e-6
        c_Tright = eos_sos_T(rho+drho,T)
        c_Tleft  = eos_sos_T(rho-drho,T)
        dc_Tdrho = (c_Tright - c_Tleft)/(np.float64(2.)*drho)
        return dc_Tdrho

def dp_TdT(rho,T):
        dT        = 1e-6
        p_Tright  = eos_p_T(rho,T+dT)
        p_Tleft   = eos_p_T(rho,T-dT)
        dp_TdT    = (p_Tright - p_Tleft)/(np.float64(2.)*dT)
        return dp_TdT

def dc_TdT(rho,T):
        dT       = 1e-6
        c_Tright = eos_sos_T(rho,T+dT)
        c_Tleft  = eos_sos_T(rho,T-dT)
        dc_TdT   = (c_Tright - c_Tleft)/(np.float64(2.)*dT)
        return dc_TdT

def Gamma(rho,p):
        one      = np.float64(1.)
        T        = eos_T(rho,p)
        dcdrho   = dc_Tdrho(rho,T)
        dpdT     = dp_TdT(rho,T)
        dcdT     = dc_TdT(rho,T)
        dcdrho_s = dcdrho + val.Zc*val.delta*T/rho**2*dpdT*dcdT
        c        = eos_sos(rho,p)
        big_gam  = one + rho/c*dcdrho_s
        return big_gam
