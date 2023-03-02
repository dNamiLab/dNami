# -----------------------------------------------------------------------------
#
# Functions used to interpolate the baseflow and bathymetry to the desired 
# grid size 
#
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------

import os
from scipy.interpolate import RegularGridInterpolator as RGI
import numpy as np
import sys

# -- Interpolate bathy
def interp_fields(nlon, nlat,fbathy='inputs/bathy.bin'):

    # =========================================================== Bathymetry 
    fname = './bathy/bathy_{}_{}.bin'.format(nlon,nlat)

    if not os.path.isfile(fname):

        with open(fbathy, 'rb') as f:
            bathy = np.load(f)


        # -- Source
        nla, nlo = bathy.shape
        lo = np.linspace(-180., 180., nlo)
        la = np.linspace(-90, 90, nla)

        fB = RGI((lo, la), bathy.T, method='linear')

        # -- Target 
        dlot = 360./nlon/2.
        dlat = 180./nlat/2.

        lot = np.linspace(-180.+dlot, 180.-dlot, nlon)
        lat = np.linspace(-90. +dlat, 90  -dlat, nlat)

        LO, LA = np.meshgrid(lot, lat, indexing='ij')

        bathy_target = fB((LO, LA)).T

        # -- Write out in correct format
        headsize = 7; head = np.empty(headsize,dtype='float64')
        head[:] = [headsize,0,0,nlat,nlon,1,1] # implicit type conversion

        with open(fname,"wb") as fh:
                np.concatenate((head,np.reshape(bathy_target,(nlat*nlon)))).tofile(fh)
        fh.closed

        print(' > bathymetry done.')

    else:
        print(' > File already exists ...')

    # =========================================================== THERMO 
    fname1 = './thermo/press_{}_{}.bin'.format(nlon,nlat)
    fname2 = './thermo/dens_{}_{}.bin'.format(nlon,nlat)

    if (not os.path.isfile(fname1)) or  (not os.path.isfile(fname2) ):

        with open('inputs/thermo.bin', 'rb') as f:
            press = np.flipud(np.load(f).T)
            dens  = np.flipud(np.load(f).T)

        # -- Source
        nla, nlo = press.shape
        lo = np.linspace(-180., 180., nlo)
        la = np.linspace(-90, 90, nla)

        fP = RGI((lo, la), press.T, method='linear')
        fD = RGI((lo, la), dens.T, method='linear')

        # -- Target 
        dlot = 360./nlon/2.
        dlat = 180./nlat/2.

        lot = np.linspace(-180.+dlot, 180.-dlot, nlon)
        lat = np.linspace(-90. +dlat, 90  -dlat, nlat)

        LO, LA = np.meshgrid(lot, lat, indexing='ij')

        press_target = fP((LO, LA)).T
        dens_target = fD((LO, LA)).T

        # -- Write out in correct format
        headsize = 7; head = np.empty(headsize,dtype='float64')
        head[:] = [headsize,0,0,nlat,nlon,1,1] # implicit type conversion

        with open(fname1,"wb") as fh:
                np.concatenate((head,np.reshape(press_target,(nlat*nlon)))).tofile(fh)
        fh.closed

        with open(fname2,"wb") as fh:
                np.concatenate((head,np.reshape(dens_target,(nlat*nlon)))).tofile(fh)
        fh.closed

        print(' > thermo done.')

    else:
        print(' > File already exists ...')

    # =========================================================== VELOCITY 
    fname1 = './velocity/utheta_{}_{}.bin'.format(nlon,nlat)
    fname2 = './velocity/uphi_{}_{}.bin'.format(nlon,nlat)

    if (not os.path.isfile(fname1)) or  (not os.path.isfile(fname2) ):

        with open('inputs/velocity.bin', 'rb') as f:
            uphi = np.flipud(np.load(f).T)
            utheta = np.flipud(np.load(f).T)

        # -- Source
        nla, nlo = utheta.shape
        lo = np.linspace(-180., 180., nlo)
        la = np.linspace(-90, 90, nla)

        fUT = RGI((lo, la), utheta.T, method='linear')
        fUP = RGI((lo, la), uphi.T, method='linear')

        # -- Target 
        dlot = 360./nlon/2.
        dlat = 180./nlat/2.

        lot = np.linspace(-180.+dlot, 180.-dlot, nlon)
        lat = np.linspace(-90. +dlat, 90  -dlat, nlat)

        LO, LA = np.meshgrid(lot, lat, indexing='ij')

        UT_target = fUT((LO, LA)).T
        UP_target = fUP((LO, LA)).T

        # -- Write out in correct format
        headsize = 7; head = np.empty(headsize,dtype='float64')
        head[:] = [headsize,0,0,nlat,nlon,1,1] # implicit type conversion

        with open(fname1,"wb") as fh:
                np.concatenate((head,np.reshape(UT_target,(nlat*nlon)))).tofile(fh)
        fh.closed

        with open(fname2,"wb") as fh:
                np.concatenate((head,np.reshape(UP_target,(nlat*nlon)))).tofile(fh)
        fh.closed

        print(' > velocity done.')

    else:
        print(' > File already exists ...')

    return 


if __name__ == "__main__":


    nlat, nlon = int(sys.argv[1]), int(sys.argv[2])

    # Make sure the output folders exits:
    if not os.path.isdir('bathy/'): os.mkdir('bathy')
    if not os.path.isdir('thermo/'): os.mkdir('thermo')
    if not os.path.isdir('velocity/'): os.mkdir('velocity')

    # Run the function
    interp_fields(nlon,nlat)



