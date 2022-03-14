# =============================================================================
# 2D euler equations in cartesian coordinates / ideal gas formulation 
# =============================================================================

# number of dimensions
dim = 2

# include transverse terms or not
iTransx = False 
iTransy = False  

# coefficients ////////////////////////////////////////////////////////////////

coefficients = {
        'delta' : 1,
        }

# unknowns to march in time ///////////////////////////////////////////////////
varname      = {'rho' : 1,
                'u'   : 2,
                'v'   : 3,
                'et'  : 4, 
                }

varsolved = ['rho','u','v','et']

# of which
consvar      = [2,3,4] # are conservative variables

# derived local variables /////////////////////////////////////////////////////

varloc = { 'e' : ' (et - 0.5_wp*(u*u + v*v) ) ',
           'p' : 'delta*rho* ( e )',
           'c' : ' ( ( 1.0_wp + delta ) * p / rho  )**0.5_wp ',
           'dedp' : '1.0_wp/delta/rho', #AT FIXED RHO
           'cp_o_av' : ' rho * (c)*(c)*(dedp)',
           }

# -- Same stored var for both

varstored    = {
         }

# names to give to the constructor ////////////////////////////////////////////
# .. for comments in the Fortran file
rhsname = {'rho' : 'd(rho)/dt',
           'u'   : 'd(rho u)/dt',
           'v'   : 'd(rho v)/dt',
           'et'   : 'd(rho et)/dt'}

# .. name tags to use for intermediate variables created by the constructor

vnamesrc_divF = {'rho'  : 'FluR',
                  'u'    : 'FluX',
                  'v'    : 'FluY',
                  'et'   : 'FluE'}

# RHS terms ///////////////////////////////////////////////////////////////////

divFx = {  'rho' : ' [rho*u                ]_1x ', 
           'u'   : ' [rho*u*u +       p    ]_1x ', 
           'v'   : ' [rho*u*v              ]_1x ', 
           'et'  : ' [u *(rho*et +       p)]_1x ' }
 
divFy = {  'rho' : ' [rho*v                ]_1y ', 
           'u'   : ' [rho*v*u              ]_1y ', 
           'v'   : ' [rho*v*v +       p    ]_1y ', 
           'et'  : ' [v *(rho*et +       p)]_1y ' }

divF = { 
        }

Ti = divFy.copy()
Tj = divFx.copy()

for key in divFx.keys():
    divF[key] = divFx[key] + ' + ' + divFy[key]

# BOUNDARY CONDITIONS ////////////////////////////////////////////////////////

locname_bc  = {'rho' : 'bc_rho',
               'u'   : 'bc_u',
               'v'   : 'bc_v',
               'et'  : 'bc_et  '}                

# ====================================================== I ======================================================================

un = ' u ' 
ut = ' v ' 
dpdn  = '  [p ]_1x  '
dundn = ' ['+un+']_1x '
dutdn = ' ['+ut+']_1x '
drdn  = '  [rho]_1x  '

# -- ASSEMBLE COEFFICIENTS

# ------------ i1
Lwavei1 = {
        'rho' : ' 0.5_wp * min( '+un+'  - (c), 0.0_wp) * ( '+dpdn+'     -          rho * (c)   * '+dundn+'  )  ' ,  
        'u'   : ' min('+un+',0.0_wp) * ( '+drdn+'                       -       1.0_wp/(c)/(c) * '+dpdn+'   )  ' ,  
        'v'   : ' min('+un+',0.0_wp) * '+dutdn+'                                                               ' ,
        'et'  : ' 0.5_wp * min( '+un+'  + (c), 0.0_wp ) * ( '+dpdn+'    +          rho * (c)   * '+dundn+'  )  ' ,  
        }

dcoefi1 = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Lwavei1['rho']+') + ('+Lwavei1['et']+' ) ) +  ('+Lwavei1['u']+') ',  
        'u'   : ' ( ('+Lwavei1['et']+') -  ('+Lwavei1['rho']+' ) )/ (c) / rho ' , 
        'v'   : ' ('+Lwavei1['v']+') ',
        'et'  : '( ('+Lwavei1['rho']+') + ('+Lwavei1['et']+' ) )' 
        }

src_phybc_wave_i1 = {
        'rho' :' ( '+dcoefi1['rho']+'  ) ',
        'u'   :' (u * ('+dcoefi1['rho']+')  + rho * ('+dcoefi1['u']+')  ) ',
        'v'   :' (v * ('+dcoefi1['rho']+') + rho * ('+dcoefi1['v']+')) ', 
        'et'  :' (et +       (p)/rho - (cp_o_av) ) * ('+dcoefi1['rho']+')  + rho * u * ('+dcoefi1['u']+') + rho* v *('+dcoefi1['v']+')   + (cp_o_av) * ('+dcoefi1['et']+')/(c)/(c) ' 
        }

# -- Save for corner
src_corner_i1 = src_phybc_wave_i1.copy()

# ------------- imax
Lwaveimax = {
        'rho' : ' 0.5_wp * max( '+un+'  - (c),0.0_wp ) * ( '+dpdn+'     -          rho * (c)   * '+dundn+'  )  ' ,  
        'u'   : ' max('+un+',0.0_wp) * ( '+drdn+'                         -       1.0_wp/(c)/(c) * '+dpdn+'   )  ' ,  
        'v'   : ' max('+un+',0.0_wp) * '+dutdn+'                                                                 ' ,
        'et'  : ' 0.5_wp * max( '+un+'  + (c),0.0_wp ) * ( '+dpdn+'     +          rho * (c)   * '+dundn+'  )  ' ,  
        }

dcoefimax = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Lwaveimax['rho']+') + ('+Lwaveimax['et']+' ) ) +  ('+Lwaveimax['u']+') ',  
        'u'   : ' ( ('+Lwaveimax['et']+') -  ('+Lwaveimax['rho']+' ) )/ (c) / rho ' , 
        'v'   : ' ('+Lwaveimax['v']+') ',
        'et'  : '( ('+Lwaveimax['rho']+') + ('+Lwaveimax['et']+' ) )' 
        }

src_phybc_wave_imax = {
        'rho' :' ( '+dcoefimax['rho']+'  ) ',
        'u'   :' ( u * ('+dcoefimax['rho']+')  + rho * ('+dcoefimax['u']+')  ) ',
        'v'   :' (v * ('+dcoefimax['rho']+') + rho * ('+dcoefimax['v']+')) ', 
        'et'  :' (et +       (p)/rho - (cp_o_av) ) * ('+dcoefimax['rho']+')  + rho * u * ('+dcoefimax['u']+') + rho* v  *('+dcoefimax['v']+')   + (cp_o_av) * ('+dcoefimax['et']+')/(c)/(c) ' 
        }

# -- Save for corner
src_corner_imax= src_phybc_wave_imax.copy()

### ---------- ADD TRANSVERSE TERMS

if iTransx:
    for key in Ti.keys():
        src_phybc_wave_i1[key]   = src_phybc_wave_i1[key]   + ' + ' + Ti[key] 
        src_phybc_wave_imax[key] = src_phybc_wave_imax[key] + ' + ' + Ti[key] 


# ====================================================== J  ======================================================================

un = ' v ' 
ut = ' u ' 
dpdn  = ' [p ]_1y  '
dundn = ' ['+un+'     ]_1y '
dutdn = ' ['+ut+'     ]_1y '
drdn  = ' [rho        ]_1y  '


# ------------- j1

Lwavej1 = {
        'rho' : ' 0.5_wp * min( '+un+'  - (c),0.0_wp ) * ( '+dpdn+'     -          rho * (c)   * '+dundn+'  )  ' ,  
        'u'   : ' min('+un+',0.0_wp) * ( '+drdn+'                         -       1.0_wp/(c)/(c) * '+dpdn+'   )  ' ,  
        'v'   : ' min('+un+',0.0_wp) * '+dutdn+'                                                                 ' ,
        'et'  : ' 0.5_wp * min( '+un+'  + (c),0.0_wp ) * ( '+dpdn+'     +          rho * (c)   * '+dundn+'  )  ' ,  
        }

dcoefj1 = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Lwavej1['rho']+') + ('+Lwavej1['et']+' ) ) +  ('+Lwavej1['u']+') ',  
        'u'   : ' ('+Lwavej1['v']+') ',
        'v'   : ' ( ('+Lwavej1['et']+') -  ('+Lwavej1['rho']+' ) )/ (c) / rho ' , 
        'et'  : '( ('+Lwavej1['rho']+') + ('+Lwavej1['et']+' ) )' 
        }

src_phybc_wave_j1 = {
        'rho' :' ( '+dcoefj1['rho']+'  ) ',
        'u'   :' ( u * ('+dcoefj1['rho']+')  + rho * ('+dcoefj1['u']+')  ) ',
        'v'   :' (v * ('+dcoefj1['rho']+') + rho * ('+dcoefj1['v']+')) ', 
        'et'  :' (et +       (p)/rho - (cp_o_av) ) * ('+dcoefj1['rho']+')  + rho * u * ('+dcoefj1['u']+') + rho* v  *('+dcoefj1['v']+')   + (cp_o_av) * ('+dcoefj1['et']+')/(c)/(c) ' 
        }


# -- Save for cornrer
src_corner_j1   = src_phybc_wave_j1.copy()

# ------------- jmax

Lwavejmax = {
        'rho' : ' 0.5_wp * max( '+un+'  - (c),0.0_wp ) * ( '+dpdn+'     -          rho * (c)   * '+dundn+'  )  ' ,  
        'u'   : ' max('+un+',0.0_wp) * ( '+drdn+'                         -       1.0_wp/(c)/(c) * '+dpdn+'   )  ' ,  
        'v'   : ' max('+un+',0.0_wp) * '+dutdn+'                                                                 ' ,
        'et'  : ' 0.5_wp * max( '+un+'  + (c),0.0_wp ) * ( '+dpdn+'     +          rho * (c)   * '+dundn+'  )  ' ,  
        }

dcoefjmax = {
        'rho' : '1.0_wp/(c)/(c) * ( ('+Lwavejmax['rho']+') + ('+Lwavejmax['et']+' ) ) +  ('+Lwavejmax['u']+') ',  
        'u'   : ' ('+Lwavejmax['v']+') ',
        'v'   : ' ( ('+Lwavejmax['et']+') -  ('+Lwavejmax['rho']+' ) )/ (c) / rho ' , 
        'et'  : '( ('+Lwavejmax['rho']+') + ('+Lwavejmax['et']+' ) )' 
        }

src_phybc_wave_jmax = {
        'rho' :' ( '+dcoefjmax['rho']+'  ) ',
        'u'   :' ( u * ('+dcoefjmax['rho']+')  + rho * ('+dcoefjmax['u']+')  ) ',
        'v'   :' (v * ('+dcoefjmax['rho']+') + rho * ('+dcoefjmax['v']+')) ', 
        'et'  :' (et +       (p)/rho - (cp_o_av) ) * ('+dcoefjmax['rho']+')  + rho * u * ('+dcoefjmax['u']+') + rho* v  *('+dcoefjmax['v']+')   + (cp_o_av) * ('+dcoefjmax['et']+')/(c)/(c) ' 
        }

# -- Save for cornrer
src_corner_jmax = src_phybc_wave_jmax.copy()

### ---------- ADD TRANSVERSE TERMS
if iTransy:
    for key in Tj.keys():
        src_phybc_wave_j1[key]   = src_phybc_wave_j1[key]   + ' + ' + Tj[key] 
        src_phybc_wave_jmax[key] = src_phybc_wave_jmax[key] + ' + ' + Tj[key] 

# ====================================================== CORNERS  ======================================================================

src_phybc_wave_i1j1     = {} 
src_phybc_wave_i1jmax   = {} 
src_phybc_wave_imaxj1   = {} 
src_phybc_wave_imaxjmax = {} 

for key in Lwavejmax.keys():
    src_phybc_wave_i1j1[key]     = src_corner_i1[key]   + ' + ' + src_corner_j1[key]
    src_phybc_wave_i1jmax[key]   = src_corner_i1[key]   + ' + ' + src_corner_jmax[key]
    src_phybc_wave_imaxj1[key]   = src_corner_imax[key] + ' + ' + src_corner_j1[key]
    src_phybc_wave_imaxjmax[key] = src_corner_imax[key] + ' + ' + src_corner_jmax[key]

