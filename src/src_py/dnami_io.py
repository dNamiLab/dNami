from dnami import np, sys

# =============================================================================
# MESH
# =============================================================================

def write_grid(tree):
        wp   = tree['misc']['working precision']
        grid = tree['grid']
        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        Lx  , Ly  , Lz   = grid['geom']['Lx']  , grid['geom']['Ly']  , grid['geom']['Lz']
        x   , y   , z    = grid['geom']['x']   , grid['geom']['y']   , grid['geom']['z']

        # Check for any non-periodic directions and extend them by 2*hlo
        bc = tree['bc']['allbc']
        hlo = tree['num']['hlo']
        if ('i1' in bc) or ('imax' in bc):
                nxgb += 2*hlo
        elif ('j1' in bc) or ('jmax' in bc):
                nygb += 2*hlo
        elif ('k1' in bc) or ('kmax' in bc):
                nzgb += 2*hlo

        # write grid to file
        if tree['mpi']['dMpi'].ioproc:
                head = np.empty(6,dtype=wp)
                head[:] = [nxgb,nygb,nzgb,Lx,Ly,Lz]
                with open('./out/axes.bin',"wb") as fh:
                        np.concatenate((head,x,y,z)).tofile(fh)
                        fh.closed

# =============================================================================
# MESSAGES
# =============================================================================

def hello_world(tree):
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']

        ioproc = dmpi.ioproc    
        iMpi   = dmpi.iMpi
        if iMpi : nprocs = dmpi.nprocs
        nxpr, nypr, nzpr  = tree['mpi']['split']['nxpr'], tree['mpi']['split']['nypr'], tree['mpi']['split']['nzpr']
        
        nxgb, nygb, nzgb  = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz    = dmpi.nx             , dmpi.ny             , dmpi.nz
        Lx  , Ly  , Lz    = grid['geom']['Lx']  , grid['geom']['Ly']  , grid['geom']['Lz']
        x   , y   , z     = grid['geom']['x']   , grid['geom']['y']   , grid['geom']['z']
        dx  , dy  , dz    = grid['geom']['dx']  , grid['geom']['dy']  , grid['geom']['dz']
        
        if ioproc:
                print('\033[1;40;94m'+'  /\     /\_   '+'\033[1;40;37m'+'  ||\  |  _       o'+'\033[1;40;31m'+'   /\      __         '+'\033[0m')       
                print('\033[1;40;94m'+'\/  \  _/   \  '+'\033[1;40;37m'+' _|| \ | /_\ |\/| |'+'\033[1;40;31m'+'  /  \    /  \        '+'\033[0m')
                print('\033[1;40;94m'+'     \/      \ '+'\033[1;40;37m'+'(_||  \|/   \|  | |'+'\033[1;40;31m'+' /    \__/    \__ v0.0'+'\033[0m')
                print('\033[1m'+'Geometry:'+'\033[0m')  
                print('-- grid size (nxgb,nygb,nzgb):',nxgb,nygb,nzgb)
                if iMpi:
                        print('-- running with mpi')
                        print('-- number of requested processes:',nxpr,'x',nypr,'x',nzpr,'=',nprocs)
                else:
                        print('-- running with NO mpi')
                print('-- local grid size (nx,ny,nz):',nx,ny,nz)
                print('-- xmin/max,dx:',x[0],x[-1],dx)
                if nygb > 1: print('-- ymin/max,dy:',y[0],y[-1],dy)
                if nzgb > 1: print('-- zmin/max,dz:',z[0],z[-1],dz)
                sys.stdout.flush()

# =============================================================================
# RESTARTS
# =============================================================================

def write_restart(n,t,flag,tree,fpath='./restarts/'):
        # file name switch (for live views)
        ndim = tree['eqns']['ndim']     
        bcs   = tree['bc']
        fnameshell = {} 
        if flag == 0:
                fname      = fpath + 'restart_' + str(n).zfill(8)
                if  len(bcs['allbc']) != 0:     
                        dirBC = ['i1','imax']
                        if ndim == 2:
                                dirBC = dirBC + ['j1','jmax']
                        if ndim == 3:   
                                dirBC = dirBC + ['j1','jmax','k1','kmax']
                        for dire in dirBC:
                                fnameshell[dire] = fpath + 'restartshell_' + str(n).zfill(8) + '_%s' % dire
                
        else:
                fname = './out/liv/restart.bin'
        
        # unpack useful tree data
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        # nvar = tree['eqns']['qvec']['nvars']
        nvar = len(tree['eqns']['qvec']['solved'])      
        hlo  = tree['num']['hlo']
        q    = tree['eqns']['qvec']['views']['q']

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz
        
        # header
        headsize = 7; head = np.empty(headsize,dtype=wp)
        head[:] = [headsize,n,t,nxgb,nygb,nzgb,nvar] # implicit type conversion
        
        # body
        dat = np.empty((nx,ny,nz,nvar),dtype=wp)
        if ndim == 3:
                dat = q[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,0:nvar].copy()
        elif ndim == 2:
                dat = q[hlo:hlo+nx,hlo:hlo+ny,np.newaxis,0:nvar].copy()
        else:
                dat = q[hlo:hlo+nx,np.newaxis,np.newaxis,0:nvar].copy()
        if iMpi:
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_WRONLY|MPI.MODE_CREATE)
                fh.Set_view(0,MPIWP,header)
                if ioproc: fh.Write_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1,0))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Write_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"wb") as fh:
                        np.concatenate((head,np.reshape(dat,(nx*ny*nz*nvar)))).tofile(fh)
                fh.closed

        # write shell (if necessary)    

        if fnameshell != {}     :
                position = {'i':[dmpi.ibeg,dmpi.iend],
                                'j':[dmpi.jbeg,dmpi.jend],
                                'k':[dmpi.kbeg,dmpi.kend]}

                nglb    = {'i':nxgb,'j':nygb,'k':nzgb}

                fh      = {}
                header  = {}
                subarray= {}

                for dir in dirBC:
                        wposition = {'i':dmpi.ibeg,
                                     'j':dmpi.jbeg,
                                     'k':dmpi.kbeg}

                        if ndim == 3:             
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb + 2*hlo}
                                extpos = ['i','j','k']
                        elif ndim == 2: 
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb }
                                extpos = ['i','j']
                        elif ndim ==1:
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb ,'k': nzgb }
                                extpos = ['i']                                          

                        # sizeglb = {'i':nxgb ,'j':nygb ,'k': nzgb }
                        sizeloc = {'i':nx  ,'j':ny  ,'k':nz  }
                        index   = {'i':[hlo,hlo+nx]  ,'j':[hlo,hlo+ny]  ,'k':[hlo,hlo+nz] }


                        for d in extpos:
                                if d != dir[0]:
                                        if position[d][0] == 1:
                                                index[d][0]= 0                                  
                                                sizeloc[d] = sizeloc[d] + hlo
                                        else:
                                                wposition[d] =  wposition[d] + hlo
                                        if position[d][1] == nglb[d]:
                                                index[d][1]= index[d][1] + hlo
                                                sizeloc[d] = sizeloc[d]  + hlo

                        sizeglb[dir[0]] = hlo

                        # header
                        headsize = 7; head = np.empty(headsize,dtype=wp)
                        head[:]  = [headsize,n,t,sizeglb['i'],sizeglb['j'],sizeglb['k'],nvar] # implicit type conversion        
                
                        if iMpi:
                                header[dir] = MPIWP.Create_contiguous(headsize)
                                header[dir].Commit()    
                                fh[dir] = MPI.File.Open(dmpi.comm_torus,fnameshell[dir],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                fh[dir].Set_view(0,MPIWP,header[dir])
                                if ioproc: fh[dir].Write_at(0,head)
                                header[dir].Free()      
                                fh[dir].Close()
                                        
                        if (position[dir[0]][0] == 1 and dir[1] == '1') or (position[dir[0]][1] == nglb[dir[0]] and dir[1:] == 'max' ):

                                sizeglb[dir[0]] = hlo
                                sizeloc[dir[0]] = hlo

                                if dir[1:] == 'max':
                                        index[dir[0]]   = [index[dir[0]][1],index[dir[0]][1]+hlo]       
                                else:                                   
                                        index[dir[0]]   = [0,hlo]                       

                                wposition[dir[0]] = 1

                                # body
                                dat     =     np.empty((sizeloc['i'],sizeloc['j'],sizeloc['k'],nvar),dtype=wp)  

                                if ndim == 3:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,index['j'][0]:index['j'][1]
                                                           ,index['k'][0]:index['k'][1],0:nvar].copy()
                                elif ndim == 2:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,index['j'][0]:index['j'][1]
                                                           ,np.newaxis,0:nvar].copy()   
                                else:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,np.newaxis
                                                           ,np.newaxis,0:nvar].copy()
                                if iMpi: 
                                        fh[dir]  = MPI.File.Open(dmpi.combc[dir],fnameshell[dir],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                        subarray = MPIWP.Create_subarray((  sizeglb['i'],    sizeglb['j'],    sizeglb['k'],  nvar),
                                                                                                         (  sizeloc['i'],    sizeloc['j'],    sizeloc['k'],  nvar),
                                                                                                         (wposition['i']-1,wposition['j']-1,wposition['k']-1,   0))
                                        subarray.Commit()
                                        disp = MPIWP.Get_size()*headsize
                                        fh[dir].Set_view(disp,MPIWP,subarray)
                                        fh[dir].Write_all(dat)
                                        subarray.Free()
                                        fh[dir].Close()
                                else:
                                        with open(fnameshell[dir],"wb") as fh[dir]:
                                                np.concatenate((head,np.reshape(dat,(sizeglb['i']*sizeglb['j']*sizeglb['k']*nvar)))).tofile(fh[dir])
                                        fh[dir].closed  


                
def read_restart(tree,fname='restart.bin'):
        # unpack useful tree data
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        # nvar = tree['eqns']['qvec']['nvars']
        nvar = len(tree['eqns']['qvec']['solved'])
        hlo  = tree['num']['hlo']
        ndim = tree['eqns']['ndim']
        q    = tree['eqns']['qvec']['views']['q']
        bcs  = tree['bc']

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz

        # read in restart
        headsize = 7
        if iMpi:
                head = np.empty(headsize,dtype=wp)
                dat  = np.empty((nx,ny,nz,nvar),dtype=wp)
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_RDONLY)
                fh.Set_view(0,MPIWP,header)
                fh.Read_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1,0))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Read_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"rb") as fh:
                        head = np.fromfile(fh,dtype=wp,count=headsize)
                        dat  = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nvar)
                        dat  = np.reshape(dat,(nx,ny,nz,nvar))
                fh.closed
        
        # run some checks & print some info to terminal:
        if ioproc:
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1;32m'+'RESTARTING FROM RESTART.BIN FILE'+'\033[0m')
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1m'+'From header:'+'\033[0m')
                print('n,time:',int(head[1]),',',head[2])
                print('nxgb,nygb,nzgb,nvar:',int(head[3]),',',int(head[4]),',',int(head[5]),',',int(head[6]))
                if nxgb != head[3]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nxgb does not match that of restart.bin'); sys.exit()
                if nygb != head[4]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nygb does not match that of restart.bin'); sys.exit()
                if nzgb != head[5]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nzgb does not match that of restart.bin'); sys.exit()
                if nvar != head[6]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nvar does not match that of restart.bin'); sys.exit()
                print('We are good to go!')
        
        # set iteration number and time:
        tree['num']['tint']['itn'] = int(head[1])
        tree['eqns']['time']       = head[2]
        
        # set fields:
        if ndim == 3:
                q[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,0:nvar] = dat.copy()
        elif ndim == 2:
                q[hlo:hlo+nx,hlo:hlo+ny,           0:nvar] = dat.reshape(nx,ny,nvar).copy()
        elif ndim == 1:
                q[hlo:hlo+nx,                      0:nvar] = dat.reshape(nx,nvar).copy()
        else:
                if ioproc: print('[error] (read_restart) ndim (',ndim,') should be 1, 2 or 3')
                sys.exit()


        # read shell (if necessary)     
        fnameshell = {} 
        if  len(bcs['allbc']) != 0:     

                dirBC = ['i1','imax']
                if ndim == 2:
                        dirBC = dirBC + ['j1','jmax']
                if ndim == 3:   
                        dirBC = dirBC + ['j1','jmax','k1','kmax']

                for dir in dirBC:
                        fnameshell[dir] = 'restartshell_'+dir

        if fnameshell != {}     :

                position = {'i':[dmpi.ibeg,dmpi.iend],
                                'j':[dmpi.jbeg,dmpi.jend],
                                'k':[dmpi.kbeg,dmpi.kend]}

                nglb    = {'i':nxgb,'j':nygb,'k':nzgb}

                fh      = {}
                header  = {}
                subarray= {}

                for dir in ['i1','imax','j1','jmax','k1','kmax']:

                        wposition = {'i':dmpi.ibeg,
                                     'j':dmpi.jbeg,
                                     'k':dmpi.kbeg}

                        if ndim == 3:             
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb + 2*hlo}
                                extpos = ['i','j','k']
                        elif ndim == 2: 
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb }
                                extpos = ['i','j']
                        elif ndim ==1:
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb ,'k': nzgb }
                                extpos = ['i']                                          

                        # sizeglb = {'i':nxgb ,'j':nygb ,'k': nzgb }
                        sizeloc = {'i':nx  ,'j':ny  ,'k':nz  }
                        index   = {'i':[hlo,hlo+nx]  ,'j':[hlo,hlo+ny]  ,'k':[hlo,hlo+nz] }


                        for d in extpos:
                                if d != dir[0]:
                                        if position[d][0] == 1:
                                                index[d][0]= 0                                  
                                                sizeloc[d] = sizeloc[d] + hlo
                                        else:
                                                wposition[d] =  wposition[d] + hlo
                                        if position[d][1] == nglb[d]:
                                                index[d][1]= index[d][1] + hlo
                                                sizeloc[d] = sizeloc[d]  + hlo


                        # header
                        headsize = 7; head = np.empty(headsize,dtype=wp)
                        head[:]  = [headsize,n,t,nxgb,nygb,nzgb,nvar] # implicit type conversion        
                
                        if iMpi:
                                header[dir] = MPIWP.Create_contiguous(headsize)
                                header[dir].Commit()    
                                fh[dir] = MPI.File.Open(dmpi.comm_torus,fnameshell[dir],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                fh[dir].Set_view(0,MPIWP,header[dir])
                                if ioproc: fh[dir].Read_at(0,head)
                                header[dir].Free()      
                                fh[dir].Close()
                                        

                        if (position[dir[0]][0] == 1 and dir[1] == '1') or (position[dir[0]][1] == nglb[dir[0]] and dir[1:] == 'max' ):

                                sizeglb[dir[0]] = hlo
                                sizeloc[dir[0]] = hlo

                                if position[dir[0]][0] == 1:
                                        index[dir[0]]   = [0,hlo]
                                else:
                                        index[dir[0]]   = [index[dir[0]][1],index[dir[0]][1]+hlo]                               

                                wposition[dir[0]] = 1

                                # body
                                dat     =     np.empty((sizeloc['i'],sizeloc['j'],sizeloc['k'],nvar),dtype=wp)  

                                if ndim == 3:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,index['j'][0]:index['j'][1]
                                                           ,index['k'][0]:index['k'][1],0:nvar].copy()
                                elif ndim == 2:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,index['j'][0]:index['j'][1]
                                                           ,np.newaxis,0:nvar].copy()   
                                else:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,np.newaxis
                                                           ,np.newaxis,0:nvar].copy()   

                                if iMpi:
                                        fh[dir]  = MPI.File.Open(dmpi.combc[dir],fnameshell[dir],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                        subarray = MPIWP.Create_subarray((  sizeglb['i'],    sizeglb['j'],    sizeglb['k'],  nvar),
                                                                                                         (  sizeloc['i'],    sizeloc['j'],    sizeloc['k'],  nvar),
                                                                                                         (wposition['i']-1,wposition['j']-1,wposition['k']-1,   0))
                                        subarray.Commit()
                                        disp = MPIWP.Get_size()*headsize
                                        fh[dir].Set_view(disp,MPIWP,subarray)
                                        fh[dir].Write_all(dat)
                                        subarray.Free()
                                        fh[dir].Close()
                                else:
                                        with open(fnameshell[dir],"wb") as fh[dir]:
                                                np.concatenate((head,np.reshape(dat,(sizeloc['i']*sizeloc['j']*sizeloc['k']*nvar)))).tofile(fh[dir])
                                        fh[dir].closed          
        return tree


# -- READ 1D BASEFLOW INTO 2D FIELD

def read_restart_1d_to_2d(tree,path=''):
        # unpack useful tree data
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        nvar = len(tree['eqns']['qvec']['solved'])
        hlo  = tree['num']['hlo']
        ndim = tree['eqns']['ndim']
        q    = tree['eqns']['qvec']['views']['q']
        bcs   = tree['bc']

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz

        # read in restart
        fname = path + 'restart.bin'
        headsize = 7

        if iMpi:
            head = np.empty(headsize,dtype=wp)
            dat  = np.empty((nx,1,1,nvar-1),dtype=wp)
            header = MPIWP.Create_contiguous(headsize)
            header.Commit() 
            fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_RDONLY)
            fh.Set_view(0,MPIWP,header)
            fh.Read_at(0,head)
            header.Free()   
            subarray = MPIWP.Create_subarray((nxgb,1,1,nvar-1),(nx,1,1,nvar-1),(0,0,0,0))
            subarray.Commit()
            disp = MPIWP.Get_size()*headsize
            fh.Set_view(disp,MPIWP,subarray)
            fh.Read_all(dat)
            subarray.Free()
            fh.Close()
        else:
            with open(fname,"rb") as fh:
                    head = np.fromfile(fh,dtype=wp,count=headsize)
                    dat  = np.fromfile(fh,dtype=wp,count=nx*1*1*(nvar-1))
                    dat  = np.reshape(dat,(nx,1,1,nvar-1))
            fh.closed

        # run some checks & print some info to terminal:
        if ioproc:
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1;32m'+'RESTARTING FROM RESTART.BIN FILE'+'\033[0m')
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1m'+'From header:'+'\033[0m')
                print('n,time:',int(head[1]),',',head[2])
                print('nxgb,nygb,nzgb,nvar:',int(head[3]),',',int(head[4]),',',int(head[5]),',',int(head[6]))
                if nxgb != head[3]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nxgb does not match that of restart.bin'); sys.exit()
                print('We are good to go!')
        
        # set iteration number and time:
        tree['num']['tint']['itn'] = int(head[1])
        tree['eqns']['time']       = head[2]
        
        # set fields:
        dat_temp = dat.reshape(nx,nvar-1).copy()  
        for k in [0,1,2]:
            for j in range(hlo,ny+hlo):
                if k == 2:
                    q[hlo:hlo+nx,j,k+1] = dat_temp[:,k] 
                else:
                    q[hlo:hlo+nx,j,k] = dat_temp[:,k] 

        return tree


# -- READ CUSTOM FILES 

def read_custom_baseflow(tree, path=None):
        # unpack useful tree data
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        nvar = 3 
        hlo  = tree['num']['hlo']
        ndim = tree['eqns']['ndim']
        q  = tree['eqns']['qvec']['views']['q']
        bcs  = tree['bc']

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz

        if path == None:
            if ioproc:
                print('MISSING PATH TO READ BASEFLOW')
                exit()


        # read in restart
        fname = path + 'restart.bin'
        headsize = 7
        if iMpi:
                head = np.empty(headsize,dtype=wp)
                dat  = np.empty((nx,ny,nz,nvar),dtype=wp)
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_RDONLY)
                fh.Set_view(0,MPIWP,header)
                fh.Read_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1,0))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Read_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"rb") as fh:
                        head = np.fromfile(fh,dtype=wp,count=headsize)
                        dat  = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nvar)
                        dat  = np.reshape(dat,(nx,ny,nz,nvar))
                fh.closed
        
        # run some checks & print some info to terminal:
        if ioproc:
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1;32m'+'RESTARTING FROM RESTART.BIN FILE'+'\033[0m')
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1m'+'From header:'+'\033[0m')
                print('n,time:',int(head[1]),',',head[2])
                print('nxgb,nygb,nzgb,nvar:',int(head[3]),',',int(head[4]),',',int(head[5]),',',int(head[6]))
                if nxgb != head[3]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nxgb does not match that of restart.bin'); sys.exit()
                if nygb != head[4]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nygb does not match that of restart.bin'); sys.exit()
                if nzgb != head[5]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nzgb does not match that of restart.bin'); sys.exit()
                if nvar != head[6]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nvar does not match that of restart.bin'); sys.exit()
                print('We are good to go!')
        
        # set iteration number and time:
        tree['num']['tint']['itn'] = int(head[1])
        tree['eqns']['time']       = head[2]
        
        # set fields:
        q[hlo:hlo+nx,                      0:nvar] = dat.reshape(nx,nvar).copy()

        # read shell (if necessary)     
        fnameshell = {} 
        if  len(bcs['allbc']) != 0:     

                dirBC = ['i1','imax']
                if ndim == 2:
                        dirBC = dirBC + ['j1','jmax']
                if ndim == 3:   
                        dirBC = dirBC + ['j1','jmax','k1','kmax']

                for dir in dirBC:
                        fnameshell[dir] = path + 'restartshell_'+dir

        if fnameshell != {}     :

                position = {'i':[dmpi.ibeg,dmpi.iend],
                                'j':[dmpi.jbeg,dmpi.jend],
                                'k':[dmpi.kbeg,dmpi.kend]}

                nglb    = {'i':nxgb,'j':nygb,'k':nzgb}

                fh      = {}
                header  = {}
                subarray= {}

                for dir in ['i1','imax']:

                        wposition = {'i':dmpi.ibeg,
                                     'j':dmpi.jbeg,
                                     'k':dmpi.kbeg}

                        if ndim == 3:             
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb + 2*hlo}
                                extpos = ['i','j','k']
                        elif ndim == 2: 
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb }
                                extpos = ['i','j']
                        elif ndim ==1:
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb ,'k': nzgb }
                                extpos = ['i']                                          

                        # sizeglb = {'i':nxgb ,'j':nygb ,'k': nzgb }
                        sizeloc = {'i':nx  ,'j':ny  ,'k':nz  }
                        index   = {'i':[hlo,hlo+nx]  ,'j':[hlo,hlo+ny]  ,'k':[hlo,hlo+nz] }


                        for d in extpos:
                                if d != dir[0]:
                                        if position[d][0] == 1:
                                                index[d][0]= 0                                  
                                                sizeloc[d] = sizeloc[d] + hlo
                                        else:
                                                wposition[d] =  wposition[d] + hlo
                                        if position[d][1] == nglb[d]:
                                                index[d][1]= index[d][1] + hlo
                                                sizeloc[d] = sizeloc[d]  + hlo


                        # header
                        headsize = 7; head = np.empty(headsize,dtype=wp)
                        head[:]  = [headsize,0,0.,nxgb,nygb,nzgb,nvar] # implicit type conversion       
                
                        if iMpi:
                                header[dir] = MPIWP.Create_contiguous(headsize)
                                header[dir].Commit()    
                                fh[dir] = MPI.File.Open(dmpi.comm_torus,fnameshell[dir],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                fh[dir].Set_view(0,MPIWP,header[dir])
                                #if ioproc: fh[dir].Read_at(0,head)
                                header[dir].Free()      
                                fh[dir].Close()
                                        

                        if (position[dir[0]][0] == 1 and dir[1] == '1') or (position[dir[0]][1] == nglb[dir[0]] and dir[1:] == 'max' ):

                                sizeglb[dir[0]] = hlo
                                sizeloc[dir[0]] = hlo

                                if position[dir[0]][0] == 1:
                                        index[dir[0]]   = [0,hlo]
                                else:
                                        index[dir[0]]   = [index[dir[0]][1],index[dir[0]][1]+hlo]                               

                                wposition[dir[0]] = 1

                                # body
                                dat     =     np.empty((sizeloc['i'],sizeloc['j'],sizeloc['k'],nvar),dtype=wp)  

                                if ndim == 3:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,index['j'][0]:index['j'][1]
                                                           ,index['k'][0]:index['k'][1],0:nvar].copy()
                                elif ndim == 2:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,index['j'][0]:index['j'][1]
                                                           ,np.newaxis,0:nvar].copy()   
                                else:
                                        dat     =     q[index['i'][0]:index['i'][1]
                                                           ,np.newaxis
                                                           ,np.newaxis,0:nvar].copy()   

                                if iMpi:
                                        fh[dir]  = MPI.File.Open(dmpi.combc[dir],fnameshell[dir],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                        subarray = MPIWP.Create_subarray((  sizeglb['i'],    sizeglb['j'],    sizeglb['k'],  nvar),
                                                                                                         (  sizeloc['i'],    sizeloc['j'],    sizeloc['k'],  nvar),
                                                                                                         (wposition['i']-1,wposition['j']-1,wposition['k']-1,   0))
                                        subarray.Commit()
                                        disp = MPIWP.Get_size()*headsize
                                        fh[dir].Set_view(disp,MPIWP,subarray)
                                        fh[dir].Write_all(dat)
                                        subarray.Free()
                                        fh[dir].Close()
                                else:
                                        with open(fnameshell[dir],"wb") as fh[dir]:
                                                np.concatenate((head,np.reshape(dat,(sizeloc['i']*sizeloc['j']*sizeloc['k']*nvar)))).tofile(fh[dir])
                                        fh[dir].closed          
        return tree



# -- WRITE CUSTOM ARRAY
def write_custom(fname,tree,var2write):

        # unpack useful tree data
        ndim = tree['eqns']['ndim']     
        bcs   = tree['bc']
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        hlo  = tree['num']['hlo']
        q    = tree['eqns']['qvec']['views'][var2write]

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz
        
        # header
        headsize = 4; head = np.empty(headsize,dtype=wp)
        head[:] = [headsize,nxgb,nygb,nzgb] # implicit type conversion
        
        # body
        dat = np.empty((nx,ny,nz),dtype=wp)
        if ndim == 3:
                dat = q[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz].copy()
        elif ndim == 2:
                dat = q[hlo:hlo+nx,hlo:hlo+ny,np.newaxis].copy()
        else:
                dat = q[hlo:hlo+nx,np.newaxis,np.newaxis].copy()
        if iMpi:
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_WRONLY|MPI.MODE_CREATE)
                fh.Set_view(0,MPIWP,header)
                if ioproc: fh.Write_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb),(nx,ny,nz),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Write_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"wb") as fh:
                        np.concatenate((head,np.reshape(dat,(nx*ny*nz)))).tofile(fh)
                fh.closed

def write_custom_val(fname,tree,q):

        # unpack useful tree data
        ndim = tree['eqns']['ndim']     
        bcs   = tree['bc']
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        hlo  = tree['num']['hlo']

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz
        
        # header
        headsize = 4; head = np.empty(headsize,dtype=wp)
        head[:] = [headsize,nxgb,nygb,nzgb] # implicit type conversion
        
        # body
        dat = np.empty((nx,ny,nz),dtype=wp)
        if ndim == 3:
                dat = q[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz].copy()
        elif ndim == 2:
                dat = q[hlo:hlo+nx,hlo:hlo+ny,np.newaxis].copy()
        else:
                dat = q[hlo:hlo+nx,np.newaxis,np.newaxis].copy()
        if iMpi:
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_WRONLY|MPI.MODE_CREATE)
                fh.Set_view(0,MPIWP,header)
                if ioproc: fh.Write_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb),(nx,ny,nz),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Write_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"wb") as fh:
                        np.concatenate((head,np.reshape(dat,(nx*ny*nz)))).tofile(fh)
                fh.closed

def read_custom(fname, tree, var2read):

        # unpack useful tree data
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        hlo  = tree['num']['hlo']
        ndim = tree['eqns']['ndim']
        q    = tree['eqns']['qvec']['views'][var2read]

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz

        # read in restart
        headsize = 4
        if iMpi:
                head = np.empty(headsize,dtype=wp)
                dat  = np.empty((nx,ny,nz),dtype=wp)
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_RDONLY)
                fh.Set_view(0,MPIWP,header)
                fh.Read_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb),(nx,ny,nz),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Read_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"rb") as fh:
                        head = np.fromfile(fh,dtype=wp,count=headsize)
                        dat  = np.fromfile(fh,dtype=wp,count=nx*ny*nz)
                        dat  = np.reshape(dat,(nx,ny,nz))
                fh.closed
        
        # run some checks & print some info to terminal:
        if ioproc:
                print('\033[1;32m'+'====================================='+'\033[0m')
                print('\033[1;32m'+'  READING IN CUSTOM VALUES TO FILL '+'\033[0m')
                print('\033[1;32m'+'    Current variable: ', var2read,' '+'\033[0m')
                print('\033[1;32m'+'====================================='+'\033[0m')
                print('\033[1m'+'From header:'+'\033[0m')
                print('nxgb,nygb,nzgb:',int(head[1]),',',int(head[2]),',',int(head[3]))
                if nxgb != head[1]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nxgb does not match that of restart.bin'); sys.exit()
                if nygb != head[2]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nygb does not match that of restart.bin'); sys.exit()
                if nzgb != head[3]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nzgb does not match that of restart.bin'); sys.exit()
                print('We are good to go!')
        
        # set fields:
        if ndim == 3:
                q[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz] = dat.copy()
        elif ndim == 2:
                q[hlo:hlo+nx,hlo:hlo+ny           ] = dat.reshape(nx,ny).copy()
        elif ndim == 1:
                q[hlo:hlo+nx                      ] = dat.reshape(nx).copy()
        else:
                if ioproc: print('[error] (read_restart) ndim (',ndim,') should be 1, 2 or 3')
                sys.exit()


        return tree
# =============================================================================
# Output
# =============================================================================

def write_data(field,n,t,tree,fpath='./out/',fname='output'):
        # File name switch (for live views)
        fprefix = fname
        ndim = tree['eqns']['ndim']     
        bcs  = tree['bc']
        fname      = fpath + '%s_' % fprefix + str(n).zfill(8)
        fnameshell = {} 
        if  len(bcs['allbc']) != 0:     
                dirBC = ['i1','imax']
                if ndim == 2:
                        dirBC = dirBC + ['j1','jmax']
                if ndim == 3:   
                        dirBC = dirBC + ['j1','jmax','k1','kmax']
                for dire in dirBC:
                        fnameshell[dire] = fpath + '%sshell_' % fprefix + dire +'_'+ str(n).zfill(8)
        
        
        # unpack useful tree data
        wp   = tree['misc']['working precision']
        dmpi = tree['mpi']['dMpi']
        grid = tree['grid']
        nvar = len(field)

        hlo  = tree['num']['hlo']

        ioproc = dmpi.ioproc
        iMpi   = dmpi.iMpi
        if iMpi : 
                nprocs = dmpi.nprocs
                MPIWP  = dmpi.MPIWP
                MPI    = dmpi.MPIlib

        nxgb, nygb, nzgb = grid['size']['nxgb'], grid['size']['nygb'], grid['size']['nzgb']
        nx  , ny  , nz   = dmpi.nx             , dmpi.ny             , dmpi.nz
        
        # header
        headsize = 7; head = np.empty(headsize,dtype=wp)
        head[:] = [headsize,n,t,nxgb,nygb,nzgb,nvar] # implicit type conversion
        
        # body
        if ndim == 3:
                dat = tree['eqns']['qvec']['views'][field[0]][hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,np.newaxis].copy()
                for f in field[1:]:
                        dat = np.append(dat,tree['eqns']['qvec']['views'][f][hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,np.newaxis].copy(),axis=3)        
        elif ndim == 2:
                dat = tree['eqns']['qvec']['views'][field[0]][hlo:hlo+nx,hlo:hlo+ny,np.newaxis,np.newaxis].copy()
                for f in field[1:]:
                        dat = np.append(dat,tree['eqns']['qvec']['views'][f][hlo:hlo+nx,hlo:hlo+ny,np.newaxis,np.newaxis].copy(),axis=3)        
        else:
                dat = tree['eqns']['qvec']['views'][field[0]][hlo:hlo+nx,np.newaxis,np.newaxis,np.newaxis].copy()
                for f in field[1:]:
                        dat = np.append(dat,tree['eqns']['qvec']['views'][f][hlo:hlo+nx,np.newaxis,np.newaxis,np.newaxis].copy(),axis=3)                
                
        if iMpi:
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(dmpi.comm_torus,fname,MPI.MODE_WRONLY|MPI.MODE_CREATE)
                fh.Set_view(0,MPIWP,header)
                if ioproc: fh.Write_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(dmpi.ibeg-1,dmpi.jbeg-1,dmpi.kbeg-1,0))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Write_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"wb") as fh:
                        np.concatenate((head,np.reshape(dat,(nx*ny*nz*nvar)))).tofile(fh)
                fh.closed

        # write shell (if necessary)    

        if fnameshell != {}     :
                position = {'i':[dmpi.ibeg,dmpi.iend],
                                'j':[dmpi.jbeg,dmpi.jend],
                                'k':[dmpi.kbeg,dmpi.kend]}

                nglb    = {'i':nxgb,'j':nygb,'k':nzgb}

                fh      = {}
                header  = {}
                subarray= {}

                for dire in dirBC:
                        wposition = {'i':dmpi.ibeg,
                                     'j':dmpi.jbeg,
                                     'k':dmpi.kbeg}


                        if ndim == 3:             
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb + 2*hlo}
                                extpos = ['i','j','k']
                        elif ndim == 2: 
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb + 2*hlo,'k': nzgb }
                                extpos = ['i','j']
                        elif ndim ==1:
                                sizeglb = {'i':nxgb + 2*hlo,'j':nygb ,'k': nzgb }
                                extpos = ['i']                                          

                        # sizeglb = {'i':nxgb ,'j':nygb ,'k': nzgb }
                        sizeloc = {'i':nx  ,'j':ny  ,'k':nz  }
                        index   = {'i':[hlo,hlo+nx]  ,'j':[hlo,hlo+ny]  ,'k':[hlo,hlo+nz] }


                        for d in extpos:
                                if d != dire[0]:
                                        if position[d][0] == 1:
                                                index[d][0]= 0                                  
                                                sizeloc[d] = sizeloc[d] + hlo
                                        else:
                                                wposition[d] =  wposition[d] + hlo
                                        if position[d][1] == nglb[d]:
                                                index[d][1]= index[d][1] + hlo
                                                sizeloc[d] = sizeloc[d]  + hlo                               

                        sizeglb[dire[0]] = hlo

                        # header
                        headsize = 7; head = np.empty(headsize,dtype=wp)
                        head[:]  = [headsize,n,t,sizeglb['i'],sizeglb['j'],sizeglb['k'],nvar] # implicit type conversion        
                
                        if iMpi:
                                header[dire] = MPIWP.Create_contiguous(headsize)
                                header[dire].Commit()    
                                fh[dire] = MPI.File.Open(dmpi.comm_torus,fnameshell[dire],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                fh[dire].Set_view(0,MPIWP,header[dire])
                                if ioproc: fh[dire].Write_at(0,head)
                                header[dire].Free()      
                                fh[dire].Close()
                                        
                        if (position[dire[0]][0] == 1 and dire[1] == '1') or (position[dire[0]][1] == nglb[dire[0]] and dire[1:] == 'max' ):

                                sizeglb[dire[0]] = hlo
                                sizeloc[dire[0]] = hlo

                                if dire[1:] == 'max':
                                        index[dire[0]]   = [index[dire[0]][1],index[dire[0]][1]+hlo]       
                                else:                                   
                                        index[dire[0]]   = [0,hlo]                       

                                wposition[dire[0]] = 1

                                # body
                                if ndim == 3:
                                        dat = tree['eqns']['qvec']['views'][field[0]][index['i'][0]:index['i'][1]
                                                                                             ,index['j'][0]:index['j'][1]
                                                                                             ,index['k'][0]:index['k'][1],np.newaxis].copy()
                                        for f in field[1:]:
                                                dat = np.append(dat,tree['eqns']['qvec']['views'][f][index['i'][0]:index['i'][1]
                                                                                                        ,index['j'][0]:index['j'][1]
                                                                                                        ,index['k'][0]:index['k'][1],np.newaxis].copy(),axis=3)
                                                                   
                                elif ndim == 2:
                                        dat = tree['eqns']['qvec']['views'][field[0]][index['i'][0]:index['i'][1]
                                                                                             ,index['j'][0]:index['j'][1]
                                                                                             ,np.newaxis,np.newaxis].copy()
                                        for f in field[1:]:
                                                dat = np.append(dat,tree['eqns']['qvec']['views'][f][index['i'][0]:index['i'][1]
                                                                                                        ,index['j'][0]:index['j'][1]
                                                                                                        ,np.newaxis,np.newaxis].copy(),axis=3)                                                     
                                else:
                                        dat = tree['eqns']['qvec']['views'][field[0]][index['i'][0]:index['i'][1]
                                                                                             ,np.newaxis
                                                                                             ,np.newaxis,np.newaxis].copy()
                                        for f in field[1:]:
                                                dat = np.append(dat,tree['eqns']['qvec']['views'][f][index['i'][0]:index['i'][1]
                                                                                                        ,np.newaxis
                                                                                                        ,np.newaxis,np.newaxis].copy(),axis=3)                                                             
                                if iMpi:
                                        fh[dire]  = MPI.File.Open(dmpi.combc[dire],fnameshell[dire],MPI.MODE_WRONLY|MPI.MODE_CREATE)
                                        subarray = MPIWP.Create_subarray((  sizeglb['i'],    sizeglb['j'],    sizeglb['k'],  nvar),
                                                                                                         (  sizeloc['i'],    sizeloc['j'],    sizeloc['k'],  nvar),
                                                                                                         (wposition['i']-1,wposition['j']-1,wposition['k']-1,   0))
                                        subarray.Commit()
                                        disp = MPIWP.Get_size()*headsize
                                        fh[dire].Set_view(disp,MPIWP,subarray)
                                        fh[dire].Write_all(dat)
                                        subarray.Free()
                                        fh[dire].Close()
                                else:
                                        with open(fnameshell[dire],"wb") as fh[dire]:
                                                np.concatenate((head,np.reshape(dat,(sizeglb['i']*sizeglb['j']*sizeglb['k']*nvar)))).tofile(fh[dire])
                                        fh[dire].closed 
                                        

# =============================================================================
# FOR USE IN RUNTIME INFO (IN TERMINAL)
# =============================================================================

def globalMinMax(tree,a,s):
        # get the global min/max value of array 'a' with string tag 's'
        
        minval = a.min()
        maxval = a.max()

        dmpi = tree['mpi']['dMpi']
        iMpi   = dmpi.iMpi
        ioproc = dmpi.ioproc

        if iMpi: 
                MPI    = dmpi.MPIlib
                comm   = dmpi.comm_torus
                globmin = comm.reduce(minval,op=MPI.MIN,root=0)
                globmax = comm.reduce(maxval,op=MPI.MAX,root=0)
        else:
                globmin = minval
                globmax = maxval
        
        if ioproc:
                if np.isnan(globmin) or np.isnan(globmax):
                        print('[error] NaN detected!')
                        sys.stdout.flush()
                        if iMpi:
                                comm.Abort()
                        else:
                                sys.exit()
                else:
                        print(s+'min,max:',globmin,globmax)
        return

def globalMax(tree,a,s):
        # get the global max value of array 'a' with string tag 's'
        
        maxval = a.max()

        dmpi = tree['mpi']['dMpi']
        iMpi   = dmpi.iMpi
        ioproc = dmpi.ioproc

        if iMpi: 
                MPI    = dmpi.MPIlib
                comm   = dmpi.comm_torus
                globmax = comm.reduce(maxval,op=MPI.MAX,root=0)
        else:
                globmax = maxval
        
        if ioproc: print(s+' max:',globmax)
        
        return

# FOR CONSTANT CFL RUNS

def dtMax(tree,a):
        # get the global min value of all dt's obtained at target CFL
        # this would be the dt to not exceed, i.e. dtMax
        
        minval = a.min()

        dmpi = tree['mpi']['dMpi']
        iMpi   = dmpi.iMpi
        ioproc = dmpi.ioproc

        if iMpi: 
                MPI    = dmpi.MPIlib
                comm   = dmpi.comm_torus
                globmin = comm.reduce(minval,op=MPI.MIN,root=0)
        else:
                globmin = minval
        
        max_dt = globmin

        return max_dt

def set_dt(tree,dt):
        # broadcast dt value from global min search and assign it

        numerics = tree['num']
        
        dmpi = tree['mpi']['dMpi']
        iMpi   = dmpi.iMpi
        ioproc = dmpi.ioproc

        if iMpi: 
                MPI    = dmpi.MPIlib
                comm   = dmpi.comm_torus
                dt = comm.bcast(dt,root=0)
        
        if ioproc:
                if np.isnan(dt):
                        print('[error] NaN detected!')
                        sys.stdout.flush()
                        if iMpi:
                                comm.Abort()
                        else:
                                sys.exit()
        
        numerics['tint']['tstep'] = dt
        
        return

'''
if iMpi:
        if wp == 'float64': MPIWP = MPI.DOUBLE
        if wp == 'float32': MPIWP = MPI.FLOAT

def write_restart(n,t,f,flag):
        # write bare minimum to restart a run from the instant of the function call
        if flag == 0:   
                fname = './restarts/restart_' + str(n).zfill(8)
        else:
                fname = './out/liv/restart.bin'
        # header
        headsize = 7
        head    = np.empty(headsize,dtype=wp)
        head[:] = [headsize,n,t,nxgb,nygb,nzgb,nvar] # implicit type conversion
        # body
        dat = np.empty((nx,ny,nz,nvar),dtype=wp)
        if ndim == 3:
                dat = f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,np.newaxis].copy()
        elif ndim == 2:
                dat = f[hlo:hlo+nx,hlo:hlo+ny,np.newaxis,np.newaxis].copy()
        else:
                dat = f[hlo:hlo+nx,np.newaxis,np.newaxis,np.newaxis].copy()     
        if iMpi:
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(comm_torus,fname,MPI.MODE_WRONLY|MPI.MODE_CREATE)
                fh.Set_view(0,MPIWP,header)
                if ioproc: fh.Write_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(ibeg-1,jbeg-1,kbeg-1,0))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Write_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"wb") as fh:
                        np.concatenate((head,np.reshape(dat,(nx*ny*nz*nvar)))).tofile(fh)
                fh.closed
        return

def read_restart():
        # read in restart
        fname = 'restart.bin'
        headsize = 7
        if iMpi:
                head = np.empty(headsize,dtype=wp)
                dat  = np.empty((nx,ny,nz,nvar),dtype=wp)
                header = MPIWP.Create_contiguous(headsize)
                header.Commit() 
                fh = MPI.File.Open(comm_torus,fname,MPI.MODE_RDONLY)
                fh.Set_view(0,MPIWP,header)
                fh.Read_at(0,head)
                header.Free()   
                subarray = MPIWP.Create_subarray((nxgb,nygb,nzgb,nvar),(nx,ny,nz,nvar),(ibeg-1,jbeg-1,kbeg-1,0))
                subarray.Commit()
                disp = MPIWP.Get_size()*headsize
                fh.Set_view(disp,MPIWP,subarray)
                fh.Read_all(dat)
                subarray.Free()
                fh.Close()
        else:
                with open(fname,"rb") as fh:
                        head = np.fromfile(fh,dtype=wp,count=headsize)
                        dat  = np.fromfile(fh,dtype=wp,count=nx*ny*nz*nvar)
                        dat  = np.reshape(dat,(nx,ny,nz,nvar))
                fh.closed
        # run some checks & print some info to terminal:
        if ioproc:
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1;32m'+'RESTARTING FROM RESTART.BIN FILE'+'\033[0m')
                print('\033[1;32m'+'================================'+'\033[0m')
                print('\033[1m'+'From header:'+'\033[0m')
                print('n,time:',int(head[1]),',',head[2])
                print('nxgb,nygb,nzgb,nvar:',int(head[3]),',',int(head[4]),',',int(head[5]),',',int(head[6]))
                if nxgb != head[3]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nxgb does not match that of restart.bin'); sys.exit()
                if nygb != head[4]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nygb does not match that of restart.bin'); sys.exit()
                if nzgb != head[5]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nzgb does not match that of restart.bin'); sys.exit()
                if nvar != head[6]: print('\033[1;40;31m'+'[error]'+'\033[0m'+' nvar does not match that of restart.bin'); sys.exit()
                print('We are good to go!')
        # set iteration number and time:
        n  = int(head[1])
        ti = head[2]
        # set fields:
        if ndim == 3:
                f[hlo:hlo+nx,hlo:hlo+ny,hlo:hlo+nz,np.newaxis] = dat.copy()
        elif ndim == 2:
                f[hlo:hlo+nx,hlo:hlo+ny,np.newaxis,np.newaxis] = dat.copy()
        else:
                f[hlo:hlo+nx,np.newaxis,np.newaxis,np.newaxis] = dat.copy()
        return
'''
