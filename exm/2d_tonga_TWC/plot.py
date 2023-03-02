# -----------------------------------------------------------------------------
#
# Basic plotting utility for the TWC Tonga explosion case 
#
#                                                          -sw 15-JUN-20
# -----------------------------------------------------------------------------
import numpy as np
import sys
import matplotlib.pyplot as plt
import matplotlib
import matplotlib as mpl
from matplotlib import rc
import matplotlib.colors as colors
import glob
import os
from matplotlib.colors import hsv_to_rgb
import cartopy as ctp 
import values

# -- Make sure to source env_dNami or add the post_io.py to python path
from post_io import read_restart_wshell, load_ax

# -- Function to combine pressure and eta fields in one RGB field 
def data_to_colors(etaf,pf, emax, pfmax): 
    ones = np.ones_like(etaf)
    colors = np.moveaxis(np.zeros_like([etaf] * 4), 0, -1)  

    # Clip fields to range
    etaf =  np.clip(etaf,-emax,emax)
    etaf = etaf/emax     
    sign = np.sign(etaf)
    pf = np.clip(pf,-pfmax,pfmax)/pfmax
    pf = 0.5*(pf + ones)
    V = pf + abs(etaf)
    V = np.clip(V,0.,1.)

    # Fill HSV channels
    hsv = np.zeros((np.shape(etaf)[0],np.shape(etaf)[1],3))
    hsv[:,:,0] = (120*ones-60*sign)/360 
    hsv[:,:,1] = abs(etaf)             
    hsv[:,:,2] = V 

    colors[:,:,:3] = hsv_to_rgb(hsv) 
    colors[:,:,3]  = 1.

    return colors


# -
try: 
    os.mkdir('./pics')
except:
    pass

path = './'

# Load Axis
xx,yy,nxx,nyy = load_ax('./out_2layer/axes.bin')
datlist = sorted(glob.glob(path+'restarts_2layer/restart_*'))

# Convert to lat lon
xx = xx*180/np.pi - 90.
yy = yy*180/np.pi -  180.

# Grab reference values
val = values.vals()
tref = val.tref

# Cartopy parameters
clon = -175. 
proj = ctp.crs.PlateCarree(central_longitude=clon)
crs = ctp.crs.PlateCarree()

# --------------------------------------------------
print('Creating the plot ...')

ll = 0.08 
lr = 0.08
offs = 0.95 
lw = 2
fig = plt.figure(figsize = (8,4) )
ax  = fig.add_axes([ll,lr,offs-ll,offs-lr], projection=proj)

# Run through the list of restart files (assuming the first file is the 0th restart) 
for i,datname in enumerate(datlist) :
    print('Reading the baseflow from ....', datname)

    t,q = read_restart_wshell(datname,verbose=False)

    h0   = q[:,:,0] 
    p    = q[:,:,6]

    # - Keep copy of baseflow
    if i == 0:
        href   = 1* h0 
        pref   = 1* p 

    delh   = (h0 - href)*val.lref 
    delpbar = (p  - pref)*val.pref

    hmax = 0.20 #m
    pmax = 200  #Pa

    print('Delta h   [m]:',  delh.min(), delh.max())
    print('Delta phi [Pa]:',  delpbar.min(), delpbar.max())

    field = data_to_colors(delh, delpbar, hmax, pmax)

    if i == 0:

        # -- Plot the field 
        im  = ax.imshow(field,extent=(-180,180,-90,90), origin='lower', transform=crs)

        # -- Modify gridlines
        gl = ax.gridlines(lw=.3,zorder=1.5, draw_labels=True)
        gl.top_labels  = False
        gl.right_labels= False
        ax.add_feature(ctp.feature.LAND, edgecolor='black', facecolor='beige', zorder = 3)

        # -- Limits
        ax.set(ylim=[xx[0], xx[-1]])
        ax.set(xlim=[yy[0], yy[-1]])

        # -- Add title with time and ranges
        time  = t*tref/3600 
        ttl = ax.set_title( 'Time (hh:mm): {0:02.0f}:{1:02.0f}. h scale: +/- {2:2.0f} [cm]. p scale: +/- {3:2.0f} [hPa]'.format(*divmod(time * 60, 60), hmax*100, pmax/100), pad = -2  )

        figname = 'pics/' + str(i).zfill(4) + '.png'.format(i)
        print('SAVING FIG ...', figname)

        plt.savefig(figname,format='png',dpi=200)
        print('DONE')
    else:
        figname = 'pics/' + str(i).zfill(4) + '.png'.format(i)
        
        # - Will skip existing figures in the 'pics/' folder
        if not os.path.isfile(figname):

            # -- Update title string 
            time  = t*tref/3600 
            time_ttl  ='Time (hh:mm): {0:02.0f}:{1:02.0f}. h scale: +/- {2:02.0f} [cm]. p scale: +/- {3:2.0f} [hPa]'.format(*divmod(time * 60, 60), hmax*100, pmax/100)
            ttl.set_text(time_ttl)

            # -- Update values
            im.remove()
            im  = ax.imshow(field,extent=(-180,180,-90,90), origin='lower', transform=crs)

            print('SAVING FIG ...', figname)
            plt.savefig(figname,format='png',dpi=200)
            print('DONE')

