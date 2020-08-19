import numpy as np
import os
#Function
def read_baseflow(fname = './baseflow.bin', interp=False, lxin=np.float64(1.), xg = None ):

    dtype = np.float64

    #Read field information:
    with open(fname,"rb") as f:
        f.seek(0,os.SEEK_SET)  
        nt     = np.fromfile(f, dtype = dtype, count = 1)             #01
        nxt    = np.fromfile(f, dtype = dtype, count = 1)             #02
        nvt    = np.fromfile(f, dtype = dtype, count = 1)             #03
        x      = np.fromfile(f, dtype = dtype, count = int(nxt))      #04
        q_size = nvt*nxt                                              #Calculate size 
        Q      = np.fromfile(f, dtype = dtype, count = int(q_size) )  # Q vector 
        p      = np.fromfile(f, dtype = dtype, count = int(nxt)  )    # Pressure

    Q   = np.reshape(Q,(int(nxt) ,int(nvt) ),   order = 'F')

    r   = Q[:,0]
    u   = Q[:,1]
    et  = Q[:,2]    

    #Interpolate to a finer grid if requested 
    if interp:
        from scipy.interpolate import interp1d

        #Rescale x
        x = x / lxin

        f = interp1d(x,r, kind='linear' )
        r = f(xg)

        f = interp1d(x,u, kind='linear' )
        u = f(xg)

        f = interp1d(x,et, kind='linear')
        et = f(xg)

   #     print(' Data interpolated to computational grid')

    return np.asarray([r,u,et]).T 



def read_bf_and_x(fname = './baseflow.bin'):

    dtype = np.float64

    #Read field information:
    with open(fname,"rb") as f:
        f.seek(0,os.SEEK_SET)  
        nt     = np.fromfile(f, dtype = dtype, count = 1)             #01
        nxt    = np.fromfile(f, dtype = dtype, count = 1)             #02
        nvt    = np.fromfile(f, dtype = dtype, count = 1)             #03
        x      = np.fromfile(f, dtype = dtype, count = int(nxt))      #04
        q_size = nvt*nxt                                              #Calculate size 
        Q      = np.fromfile(f, dtype = dtype, count = int(q_size) )  # Q vector 
        p      = np.fromfile(f, dtype = dtype, count = int(nxt)  )    # Pressure

    Q   = np.reshape(Q,(int(nxt) ,int(nvt) ),   order = 'F')

    print(p)
    r   = Q[:,0]
    u   = Q[:,1]
    et  = Q[:,2]

    return x, np.asarray([r,u,p]).T 
