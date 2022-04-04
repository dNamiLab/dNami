# =============================================================================
# 1D euler equations in cartesian coordiantes 
#                                                          -     sw 01-SEP-21
# =============================================================================

# number of dimensions
dim = 1

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'delta' : 1,
        'amp'   : 2,
        'omega' : 3,
        't'     : 4,
        'mub'   : 5,
        'kappa' : 6,
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                'u'   : 2,
                'et'  : 3, 
                }

varsolved = ['rho','u','et']

# of which
consvar      = [2,3] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

varloc = { 'e'       : ' (et - 0.5_wp*u*u) ',
           'p'       : 'delta*rho* ( e )',
           'T'       : 'p/rho',
           'c'       : ' ( ( 1.0_wp + delta ) * p / rho  )**0.5_wp ',
           }


# RHS terms ///////////////////////////////////////////////////////////////////

F    = {  
        'rho' : ' (rho*u         ) ', 
        'u'   : ' (rho*u*u + p   ) ', 
        'et'  : ' (u*(rho*et + p)) '
        }
 
divF = { 
        'rho' : ' [ '+F['rho']+'  ]_1x  ', 
        'u'   : ' [ '+F['u'] +'   ]_1x  ', 
        'et'  : ' [ '+F['et'] +'  ]_1x  ', 
        }

visc = {
        'u'   : ' - mub * [u]_2xx ', 
        'et'  : ' - mub * ( [ u*{u}_1x  ]_1x  ) - kappa * [T]_2xx ', 
        }

# ====================================================== I ======================================================================

# ------------ i1
Li1 = {
        'rho' : ' 0.5_wp * min( u  - (c), 0.0_wp) * (  [p]_1x      -          rho * (c)   *  [u]_1x   )  ' ,  
        'u'   : ' min(u,0.0_wp) * (  [rho]_1x                        -       1.0_wp/(c)/(c) *  [p]_1x    )  ' ,  
        'et'  : ' 0.5_wp * min( u  + (c), 0.0_wp ) * (  [p]_1x     +          rho * (c)   *  [u]_1x   )  ' ,  
        }

# -- SPECIFY INLET PROPAGATING ACOUSTIC WAVE
Li1['u']  = '  +  amp*omega*cos(omega*t) '

dcoefi1 = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Li1['rho']+') + ('+Li1['et']+' ) ) +  ('+Li1['u']+') ',  
        'u'   : ' ( ('+Li1['et']+') -  ('+Li1['rho']+' ) )/ (c) / rho ' , 
        'et'  : '( ('+Li1['rho']+') + ('+Li1['et']+' ) )' 
        }

src_phybc_wave_i1 = {
        'rho' :' ( '+dcoefi1['rho']+'  ) ',
        'u'   :' (u * ('+dcoefi1['rho']+')  + rho * ('+dcoefi1['u']+')  ) ',
        'et'  :' (et +       (p)/rho - (c*c/delta) ) * ('+dcoefi1['rho']+')  + rho * u * ('+dcoefi1['u']+')  + (c*c/delta) * ('+dcoefi1['et']+')/(c)/(c) ',  
        }

# ------------- imax
Limax = {
        'rho' : ' 0.5_wp * max( u  - (c),0.0_wp ) * (  [p]_1x      -          rho * (c)   *  [u]_1x   )  ' ,  
        'u'   : ' max(u,0.0_wp) * (  [rho]_1x                          -       1.0_wp/(c)/(c) *  [p]_1x    )  ' ,  
        'et'  : ' 0.5_wp * max( u  + (c),0.0_wp ) * (  [p]_1x      +          rho * (c)   *  [u]_1x   )  ' ,  
        }


dcoefimax = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Limax['rho']+') + ('+Limax['et']+' ) ) +  ('+Limax['u']+') ',  
        'u'   : ' ( ('+Limax['et']+') -  ('+Limax['rho']+' ) )/ (c) / rho ' , 
        'et'  : '( ('+Limax['rho']+') + ('+Limax['et']+' ) )' 
        }

src_phybc_wave_imax = {
        'rho' :' ( '+dcoefimax['rho']+'  ) ',
        'u'   :' ( u * ('+dcoefimax['rho']+')  + rho * ('+dcoefimax['u']+')  ) ',
        'et'  :' (et +       (p)/rho - (c*c/delta) ) * ('+dcoefimax['rho']+')  + rho * u * ('+dcoefimax['u']+') + (c*c/delta) * ('+dcoefimax['et']+')/(c)/(c) ',  
        }

# ===============================================================================================================================

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file

rhsname = {'rho'  : 'd(rho)/dt',
           'u'    : 'd(rho u)/dt',
           'et'   : 'd(rho et)/dt',
           }

# .. Intermediate variable names 

vnamesrc_divF = {}
vnamesrc_visc = {}

for var in varsolved:
    vnamesrc_divF[var] = 'divF_' + var
    vnamesrc_visc[var] = 'visc_' + var
               
              
