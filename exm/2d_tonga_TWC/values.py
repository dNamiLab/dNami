# -----------------------------------------------------------------------------
#
# Class containing reference values for dNami  
#
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------

import numpy as np

class vals:
    
    def __init__(self):

        # Values from ISA at 80km altitude  
        self.hatm   = 80000.0                  #[m]
        self.pbar   = 9297.70532429097         #[Pa] 
        self.rbar   = 0.12945063236647714      #[kg/m3] 
        self.ctilde = 317.10261143656646       #[m/s] 
        self.Ttilde = 250.21290672215093       #[K] 

        # Values from ERA5 Jan 15th fields 
        self.hERA5 = 50610.0                    #[m]
        self.pERA5 = 15266.845456579336         #[Pa]
        self.rERA5 = 0.21517602110822553        #[kg/m3]

        # Values chosen as reference values for dNami 
        self.g     = 9.81                       #[m/s2]
        self.lref  = 3668                       #[m]    - reference length
        self.uref  = np.sqrt(self.g*self.lref)  #[m/s]  - reference velocity 
        self.tref  = self.lref/self.uref        #[s]    - reference time
        self.rref  = 1e3                        #[kg/m3]- reference density
        self.pref  = self.rref*self.uref**2     #[Pa]   - reference pressure
        self.gamma = 1.4                        #[-]    - gamma value for air 

        # 
        self.R     = 8.314
        self.M     = 0.029
        self.Rs    = self.R/self.M

    def print_info(self):
        print(' --- Atmospheric Values ---')
        print(' Atmospheric thickness                             : {0:<8g} [m]'.format(self.hatm))
        print(' Vertical average pressure                         : {0:<8g} [Pa]'.format(self.pbar))
        print(' Vertical average density                          : {0:<8g} [kg/m3]'.format(self.rbar))
        print(' Vertical average density weighted speed of sound  : {0:<8g} [m/s]'.format(self.ctilde))
        print(' Vertical average density weighted temperature     : {0:<8g} [K]'.format(self.Ttilde))
        print(' --------------')

if __name__ == '__main__':
    val = vals()
    val.print_info()
