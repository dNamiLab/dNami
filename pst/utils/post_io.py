import numpy as np

# -- Post processing tools for dNami

# === Axis loader 
def load_ax(path,wp='float64'):
    """
    Load the axes (x in 1d, (x,y) in 2d, (x,y,z) in 3d) which are the coordinates in physical space as well as the number of points in each direction. 

    Args:
        path: The path to the axis file. Usually written out in the work directory in /path/to/wrk/out/axes.bin 
    Returns:
        If its 1D problem x,nxgb is returned
        If its 2D problem x,y,nxgb,nygb is returned
        If its 3D problem x,y,z,nxgb,nygb,nzgb is returned
    """
    with open(path,"rb") as fh:
        head = np.fromfile(fh,dtype=wp,count=6)
        nxgb, nygb, nzgb = int(head[0]), int(head[1]), int(head[2])
        Lx, Ly, Lz = head[3], head[4], head[5]; del head
        if nygb == 1: ndim = 1
        elif nzgb == 1: ndim = 2
        elif nzgb >  1: ndim = 3
        x = np.fromfile(fh,dtype=wp,count=nxgb); dx = x[1]-x[0]
        if ndim>1: y = np.fromfile(fh,dtype=wp,count=nygb); dy = y[1]-y[0]
        if ndim>2: z = np.fromfile(fh,dtype=wp,count=nzgb); dz = z[1]-z[0]
        fh.closed
    if ndim == 1:
        return x,nxgb
    elif ndim == 2:
        return x,y,nxgb,nygb
    elif ndim == 3 :
        return x,y,z,nxgb,nygb,nzgb

# === Restart loader -- core only 
def read_restart(fname,wp='float64'):
    """
    Input a restart_XXXXXXXX file and the function will return the full core of q. Can also be used to read custom outputs from the dnami_io.write_data() function. 

    Args:
        fname: The path to the restart file. Usually written out in the work directory in /path/to/wrk/restarts/restart_XXXXXXXX
    Returns:
        The timestep number n, the time t and the variables in the core of the domain q are returned
    """
    with open(fname,"rb") as fh:
        headsize = int(np.fromfile(fh,dtype=wp,count=1))
        head = np.fromfile(fh,dtype=wp,count=headsize-1)
        n = int(head[0])
        t = head[1]
        nx,ny,nz,nv = int(head[2]),int(head[3]),int(head[4]),int(head[5])
        f = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nv)
        if ny == 1: ndim = 1
        elif nz == 1: ndim = 2
        elif nz >  1: ndim = 3
        if ndim == 3:
            f = np.reshape(f,(nx,ny,nz,nv))
        elif ndim == 2:
            f = np.reshape(f,(nx,ny   ,nv))
        else:
            f = np.reshape(f,(nx      ,nv))
    fh.closed
    return n,t,f

# === Restart loader -- core and shells 
def read_restart_wshell(fname,verbose=False,wp='float64'):
    """
    Input a restart_XXXXXXXX file and the function will return the full q including the shell information. Can also be used to read custom outputs from the dnami_io.write_data() function. 

    Args:
        fname: The path to the restart file. Usually written out in the work directory in /path/to/wrk/restarts/restart_XXXXXXXXX. The shell file name will be automatically detected. 
        verbose: print additional information
    Returns:
        The timestep number n, the time t and the variables in the full domain, including the shells, q are returned
    """

    if verbose:
        print('Input path:', fname)

    # -- Split path to insert 'shell' and '[nit]' 
    psplit  = fname.split('/') 
    prepath = '/'.join(psplit[:-1])
    fnamesh = psplit[-1].split('_')
    fnamesh[0] = prepath + '/' + ''.join(fnamesh[:-1])
    fnamesh[1] = fnamesh[-1]

    try:
        n, t, qi1    = read_restart(fname=fnamesh[0] +'shell_' + fnamesh[1] + '_i1'  )
        n, t, qimax  = read_restart(fname=fnamesh[0] +'shell_' + fnamesh[1] + '_imax')
        hloi = True
        if verbose:
            print(' > Read ishell ..')
    except Exception as e:
        #If no ishell, assume periodic in i
        hloi = False 

    try:
        n, t, qj1    = read_restart(fname=fnamesh[0] +'shell_' + fnamesh[1] + '_j1'  )
        n, t, qjmax  = read_restart(fname=fnamesh[0] +'shell_' + fnamesh[1] + '_jmax'  )
        hloj = True
        if verbose:
            print(' > Read jshell ..')
    except Exception as e:
        #If no jshell, assume periodic in j
        hloj = False 

    try:
        n, t, qk1    = read_restart(fname=fnamesh[0] +'shell_' + fnamesh[1] + '_k1'  )
        n, t, qkmax  = read_restart(fname=fnamesh[0] +'shell_' + fnamesh[1] + '_kmax'  )
        hlok = True
        if verbose:
            print(' > Read kshell ..')
    except Exception as e:
        hlok = False 

    n, t, qcore  = read_restart(fname=fname)

    # -- Get the dimension from the shape of qcore
    ndim = len(qcore.shape) - 1

    # -- Get hlo size
    if hloi:   hlo = qi1.shape[0]
    elif hloj: hlo = qj1.shape[1]
    elif hlok: hlo = qk1.shape[2]
    else: return t,qcore

    if verbose:
        print('Dimensions: ', ndim)

    sl_j  = None; sl_k  = None
    corej = None; corek = None

    #  -- Get core shape
    if ndim ==1:
        nx,nvar = qcore.shape
    elif ndim ==2:
        nx,ny,nvar = qcore.shape
    else :
        nx,ny,nz,nvar = qcore.shape

    #  -- Create all the necessary slices 
    if hloi:
        nxgb = nx + 2*hlo
        sl_i1   = np.s_[0:hlo]
        sl_imax = np.s_[nx+hlo:nx+2*hlo]
        corei   = np.s_[hlo:nx+hlo]
    else:
        nxgb = 1*nx
        corei = np.s_[0:nx]
    sl_i = np.s_[0:nxgb]
    if sl_i == corei: sl_ic = np.s_[hlo:nx+hlo]
    else:             sl_ic = np.s_[0:nx+2*hlo]

    if ndim >1:
        if hloj:
            nygb    = ny + 2*hlo
            sl_j1   = np.s_[0:hlo]
            sl_jmax = np.s_[ny+hlo:ny+2*hlo]
            corej   = np.s_[hlo:ny+hlo]
        else:
            nygb  = 1*ny
            corej = np.s_[0:ny]
        sl_j = np.s_[0:nygb]
        if sl_j == corej: sl_jc = np.s_[hlo:ny+hlo]
        else:             sl_jc = np.s_[0:ny+2*hlo]

    if ndim>2:
        if hlok:
            nzgb    = nz + 2*hlo
            sl_k1   = np.s_[0:hlo]
            sl_kmax = np.s_[nz+hlo:nz+2*hlo]
            corek   = np.s_[hlo:nz+hlo]
        else:
            nzgb = 1*nz
            corek = np.s_[0:nz]
        sl_k = np.s_[0:nzgb]
        sl_kc = np.s_[hlo:nz+hlo]
        if sl_k == corek: sl_kc = np.s_[hlo:nz+hlo]
        else:             sl_kc = np.s_[0:nz+2*hlo]

    # -- Shape:
    if ndim ==1:    shape = (nxgb,nvar)
    elif ndim ==2:  shape = (nxgb,nygb,nvar)
    else:           shape = (nxgb,nygb,nzgb,nvar)

    # Allocate memory
    q = np.empty(shape=shape,dtype = wp)

    if verbose:
        print('q with shells shape:' , shape)

    # Prepare some slices
    sl_var = np.s_[0:nvar]
    sl_hlo = np.s_[0:hlo]

    # x halos:
    if hloi:
        if verbose:
            print('Filling i shells')
        if ndim ==1:
            q[sl_i1  , sl_var] = qi1[sl_hlo  , sl_var]
            q[sl_imax, sl_var] = qimax[sl_hlo, sl_var]
        elif ndim ==2: 
            q[sl_i1  ,sl_j, sl_var] = qi1  [sl_hlo ,sl_jc, sl_var]
            q[sl_imax,sl_j, sl_var] = qimax[sl_hlo ,sl_jc, sl_var]
        else:
            q[sl_i1  ,sl_j, sl_k, sl_var] = qi1  [sl_hlo,sl_jc, sl_kc, sl_var]
            q[sl_imax,sl_j, sl_k, sl_var] = qimax[sl_hlo,sl_jc, sl_kc, sl_var]

    # y halos:
    if hloj:
        if verbose:
            print('Filling j shells')
        if ndim ==2: 
            q[sl_i  ,sl_j1  , sl_var] = qj1  [sl_ic, sl_hlo, sl_var]
            q[sl_i  ,sl_jmax, sl_var] = qjmax[sl_ic, sl_hlo, sl_var]
        else:
            q[sl_i  ,sl_j1  , sl_k, sl_var] = qj1  [sl_ic,sl_hlo,sl_kc,sl_var]
            q[sl_i  ,sl_jmax, sl_k, sl_var] = qjmax[sl_ic,sl_hlo,sl_kc,sl_var]

    # z halos:
    if hlok:
        if verbose:
            print('Filling k shells')
        q[sl_i  ,sl_j   , sl_k1  , sl_var] = qk1  [sl_ic,sl_jc, sl_hlo, sl_var]
        q[sl_i  ,sl_j   , sl_kmax, sl_var] = qkmax[sl_ic,sl_jc, sl_hlo, sl_var]

    # CORE:
    if ndim ==1:
        q[corei, sl_var] = qcore
    elif ndim ==2:
        q[corei, corej, sl_var] = qcore
    else:
        q[corei, corej, corek, sl_var] = qcore

    return t,q
