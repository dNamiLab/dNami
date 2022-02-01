# =============================================================================
# 1D euler equations in cartesian coordinates ig
# =============================================================================

# number of dimensions
dim = 1

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'delta' : 1,
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

varloc = { 'e' : ' (et - 0.5_wp*u*u) ',                        #internal energy
           'p' : 'delta*rho* ( e )',                           #pressure equation of state
           'c' : '( ( 1.0_wp + delta ) * p / rho  )**0.5_wp ', #isentropic speed of sound
           }

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho'  : 'd(rho)/dt',
           'u'    : 'd(rho u)/dt',
           'et'   : 'd(rho et)/dt',
           }

# .. name tags to use for intermediate variables created by the constructor
vnamesrc_divF = {'rho'  : 'FluRx',
                 'u'    : 'FluXx',
                 'et'   : 'FluEx'}

# RHS terms ///////////////////////////////////////////////////////////////////
# .. Flux vector
divF    = {  
        'rho' : ' [ rho*u          ]_1x ', 
        'u'   : ' [ rho*u*u + p    ]_1x ', 
        'et'  : ' [ u*(rho*et + p) ]_1x ', 
        }
 

# =================================================Boundary conditions==============================================================

# -- ASSEMBLE COEFFICIENTS

# ------------ i1
Lwavei1 = {
        'rho' : ' 0.5_wp * min( u  - (c), 0.0_wp) * ( [p]_1x      -          rho * (c)   * [u]_1x   )  ' ,  
        'u'   : ' min(u,0.0_wp) * ( [rho]_1x                     -       1.0_wp/(c)/(c) * [p]_1x    )  ' ,  
        'et'  : ' 0.5_wp * min( u  + (c), 0.0_wp ) * ( [p]_1x     +          rho * (c)   * [u]_1x   )  ' ,  
        }

dcoefi1 = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Lwavei1['rho']+') + ('+Lwavei1['et']+' ) ) +  ('+Lwavei1['u']+') ',  
        'u'   : ' ( ('+Lwavei1['et']+') -  ('+Lwavei1['rho']+' ) )/ (c) / rho ' , 
        'et'  : '( ('+Lwavei1['rho']+') + ('+Lwavei1['et']+' ) )' 
        }

src_phybc_wave_i1 = {
        'rho' :' ( '+dcoefi1['rho']+'  ) ',
        'u'   :' (u * ('+dcoefi1['rho']+')  + rho * ('+dcoefi1['u']+')  ) ',
        'et'  :' (et +       (p)/rho - (c*c/delta) ) * ('+dcoefi1['rho']+')  + rho * u * ('+dcoefi1['u']+')  + (c*c/delta) * ('+dcoefi1['et']+')/(c)/(c) ',  
        }

# ------------- imax
Lwaveimax = {
        'rho' : ' 0.5_wp * max( u  - (c),0.0_wp ) * ( [p]_1x      -          rho * (c)   * [u]_1x   )  ' ,  
        'u'   : ' max(u,0.0_wp) * ( [rho]_1x                      -       1.0_wp/(c)/(c) * [p]_1x   )  ' ,  
        'et'  : ' 0.5_wp * max( u  + (c),0.0_wp ) * ( [p]_1x      +          rho * (c)   * [u]_1x   )  ' ,  
        }


dcoefimax = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Lwaveimax['rho']+') + ('+Lwaveimax['et']+' ) ) +  ('+Lwaveimax['u']+') ',  
        'u'   : ' ( ('+Lwaveimax['et']+') -  ('+Lwaveimax['rho']+' ) )/ (c) / rho ' , 
        'et'  : '( ('+Lwaveimax['rho']+') + ('+Lwaveimax['et']+' ) )' 
        }

src_phybc_wave_imax = {
        'rho' :' ( '+dcoefimax['rho']+'  ) ',
        'u'   :' ( u * ('+dcoefimax['rho']+')  + rho * ('+dcoefimax['u']+')  ) ',
        'et'  :' (et +       (p)/rho - (c*c/delta) ) * ('+dcoefimax['rho']+')  + rho * u * ('+dcoefimax['u']+') + (c*c/delta) * ('+dcoefimax['et']+')/(c)/(c) ',  
        }
