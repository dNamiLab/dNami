"""
The dnami mpi module provides functions
for MPI communication.
"""

from dnami import np, sys, dnamiF
from timeit import default_timer as timer
#import os
#import time
# TODO: make swap compatible with varstored


class type_mpi:
    """
    The type_mpi object contains functions for swapping halo data among different MPI processes.

    Args:
      tree: The dnami tree data structure
    """
    # mpi class

    def __init__(self, tree):
        wp = tree['misc']['working precision']
        nxgb, nygb, nzgb = tree['grid']['size']['nxgb'], tree['grid']['size']['nygb'], tree['grid']['size']['nzgb']
        nxpr, nypr, nzpr = tree['mpi']['split']['nxpr'], tree['mpi']['split']['nypr'], tree['mpi']['split']['nzpr']
        #nvar = tree['eqns']['qvec']['nvars']
        #hlo  = tree['num']['deriv']['hlo']
        # start mpi (or not...)
        if nxpr == 1 and nypr == 1 and nzpr == 1:
            self.iMpi = False  # running with no mpi
        else:
            self.iMpi = True
        if self.iMpi:
            from mpi4py import MPI  # [note] mpi_init is done in 'import mpi'
            self.MPIlib = MPI
            # create torus
            self.comm_torus = MPI.COMM_WORLD.Create_cart(
                [nxpr, nypr, nzpr], periods=True, reorder=True)
            self.nprocs = self.comm_torus.Get_size()
            self.procid = self.comm_torus.Get_rank()
            self.coords = self.comm_torus.Get_coords(self.procid)
            # assign who is in charge of console I/Os
            self.ioproc = False
            if self.procid == 0:
                self.ioproc = True  # do not change ioproc rank value, 0 is assumed elsewhere
            # check requested number of cores makes sense
            if nxpr * nypr * nzpr != self.nprocs:
                if self.ioproc:
                    print(
                        '[error] nxpr*nypr*nzpr does not match the number of cores requested')
                sys.exit()
            # check that nxgb,nygb,nzgb are multiples of nxpr,nypr,nzpr
            if (nxgb % nxpr) == 0:
                self.nx = int(nxgb / nxpr)
            else:
                if self.ioproc:
                    print('[error] nxgb is not a multiple of nxpr')
                sys.exit()
            if (nygb % nypr) == 0:
                self.ny = int(nygb / nypr)
            else:
                if self.ioproc:
                    print('[error] nygb is not a multiple of nypr')
                sys.exit()
            if (nzgb % nzpr) == 0:
                self.nz = int(nzgb / nzpr)
            else:
                if self.ioproc:
                    print('[error] nzgb is not a multiple of nzpr')
                sys.exit()
            # create global/local index mapping
            i = 1
            for itor in range(0, nxpr):
                if self.coords[0] == itor:
                    self.ibeg = i
                    self.iend = self.ibeg + self.nx - 1
                i = i + self.nx
            i = 1
            for itor in range(0, nypr):
                if self.coords[1] == itor:
                    self.jbeg = i
                    self.jend = self.jbeg + self.ny - 1
                i = i + self.ny
            i = 1
            for itor in range(0, nzpr):
                if self.coords[2] == itor:
                    self.kbeg = i
                    self.kend = self.kbeg + self.nz - 1
                i = i + self.nz
            # create neighbour map for swapping processes
            dummy, self.neighxm = self.comm_torus.Shift(0, -1)
            dummy, self.neighxp = self.comm_torus.Shift(0, 1)
            dummy, self.neighym = self.comm_torus.Shift(1, -1)
            dummy, self.neighyp = self.comm_torus.Shift(1, 1)
            dummy, self.neighzm = self.comm_torus.Shift(2, -1)
            dummy, self.neighzp = self.comm_torus.Shift(2, 1)
            if tree['misc']['verbose']:
                self.showTorus()
            # create reference to MPI data type
            if wp == 'float64':
                self.MPIWP = MPI.DOUBLE
            if wp == 'float32':
                self.MPIWP = MPI.FLOAT

            # create boundary communicators (if needed):
            if len(tree['bc']['allbc']) != 0:

                ndim = tree['eqns']['ndim']

                edges = {'i': [1, nxgb], 'j': [1, nygb], 'k': [1, nzgb]}

                locedges = {'i': [self.ibeg, self.iend],
                            'j': [self.jbeg, self.jend],
                            'k': [self.kbeg, self.kend]}

                members = {'i1': np.empty(self.nprocs, dtype='i'),
                           'imax': np.empty(self.nprocs, dtype='i'),
                           'j1': np.empty(self.nprocs, dtype='i'),
                           'jmax': np.empty(self.nprocs, dtype='i'),
                           'k1': np.empty(self.nprocs, dtype='i'),
                           'kmax': np.empty(self.nprocs, dtype='i')}

                group = {
                    'i1': None,
                    'imax': None,
                    'j1': None,
                    'jmax': None,
                    'k1': None,
                    'kmax': None}
                commbc = {
                    'i1': None,
                    'imax': None,
                    'j1': None,
                    'jmax': None,
                    'k1': None,
                    'kmax': None}

                dirBC = ['i1', 'imax']
                if ndim == 2:
                    dirBC = dirBC + ['j1', 'jmax']
                if ndim == 3:
                    dirBC = dirBC + ['j1', 'jmax', 'k1', 'kmax']

                for dir in dirBC:

                    ibelong = 0
                    # count members
                    if dir[1] == '1':
                        if locedges[dir[0]][0] == edges[dir[0]][0]:
                            ibelong = 1
                    elif dir[1:] == 'max':
                        if locedges[dir[0]][1] == edges[dir[0]][1]:
                            ibelong = 1

                    if self.ioproc:
                        members[dir][0] = ibelong
                        for pro in range(1, self.nprocs):
                            buf = np.empty(1, dtype='i')
                            members[dir][pro] = self.comm_torus.recv(
                                source=pro, tag=pro)
                    else:
                        buf = np.empty(1, dtype='i')
                        buf[0] = ibelong
                        self.comm_torus.send(ibelong, dest=0, tag=self.procid)
                    members[dir] = self.comm_torus.bcast(members[dir], root=0)

                    gmember = [
                        pro for pro in range(
                            0, self.nprocs) if members[dir][pro] == 1]

                    group[dir] = self.comm_torus.group.Incl(gmember)
                    commbc[dir] = self.comm_torus.Create(group[dir])

                self.combc = commbc

        else:
            self.nx = nxgb
            self.ny = nygb
            self.nz = nzgb
            self.ibeg = 1
            self.iend = self.nx
            self.jbeg = 1
            self.jend = self.ny
            self.kbeg = 1
            self.kend = self.nz
            self.ioproc = True
        self.iSwap = {'i1': True, 'imax': True,
                                  'j1': True, 'jmax': True,
                                  'k1': True, 'kmax': True}

    def showTorus(self):
        # visualise torus
        if self.ioproc:
            print(
                '\033[1;104m' +
                'proc id|   with coords  | ibeg jbeg kbeg |      x-,     y-,     z-  neighbours' +
                '\033[0m')
            print(
                '\033[1;104m' +
                '       |                | iend jend kend |      x+,     y+,     z+   (proc id)' +
                '\033[0m')
            print(
                "\033[1;30;47m{0:7d}|({1:4d},{2:4d},{3:4d})| {4:4d} {5:4d} {6:4d} | {7:7d},{8:7d},{9:7d}\033[0m".format(
                    self.procid,
                    self.coords[0],
                    self.coords[1],
                    self.coords[2],
                    self.ibeg,
                    self.jbeg,
                    self.kbeg,
                    self.neighxm,
                    self.neighym,
                    self.neighzm))
            print(
                "\033[1;30;47m       |                | {0:4d} {1:4d} {2:4d} | {3:7d},{4:7d},{5:7d}\033[0m".format(
                    self.iend,
                    self.jend,
                    self.kend,
                    self.neighxp,
                    self.neighyp,
                    self.neighzp))
            for pro in range(1, self.nprocs):
                data = np.empty(16, dtype='i')
                # self.comm_torus.Recv([data,MPI.INT],source=pro,tag=pro)
                self.comm_torus.Recv(data, source=pro, tag=pro)
                if (pro % 2) == 0:
                    print("\033[1;30;47m{0:7d}|({1:4d},{2:4d},{3:4d})| {4:4d} {5:4d} {6:4d} | {7:7d},{8:7d},{9:7d}\033[0m".format(
                        data[0], data[1], data[2], data[3], data[4], data[6], data[8], data[10], data[12], data[14]))
                    print(
                        "\033[1;30;47m       |                | {0:4d} {1:4d} {2:4d} | {3:7d},{4:7d},{5:7d}\033[0m".format(
                            data[5],
                            data[7],
                            data[9],
                            data[11],
                            data[13],
                            data[15]))
                else:
                    print("\033[1;37;40m{0:7d}|({1:4d},{2:4d},{3:4d})| {4:4d} {5:4d} {6:4d} | {7:7d},{8:7d},{9:7d}\033[0m".format(
                        data[0], data[1], data[2], data[3], data[4], data[6], data[8], data[10], data[12], data[14]))
                    print(
                        "\033[1;37;40m       |                | {0:4d} {1:4d} {2:4d} | {3:7d},{4:7d},{5:7d}\033[0m".format(
                            data[5],
                            data[7],
                            data[9],
                            data[11],
                            data[13],
                            data[15]))
        else:
            data = np.empty(16, dtype='i')
            data[:] = [
                self.procid,
                self.coords[0],
                self.coords[1],
                self.coords[2],
                self.ibeg,
                self.iend,
                self.jbeg,
                self.jend,
                self.kbeg,
                self.kend,
                self.neighxm,
                self.neighxp,
                self.neighym,
                self.neighyp,
                self.neighzm,
                self.neighzp]
            # self.comm_torus.Send([data,MPI.INT],dest=0,tag=self.procid)
            self.comm_torus.Send(data, dest=0, tag=self.procid)

    # swap "conceptual" sketch with references to local python indices:
    #
    #                     negative dir swap --        ++ positive dir swap
    # PROC I-1 =============================||        ||============================= PROC I+1
    #   domain -------------->|<---halo---->||        ||<----halo--->|<-------------- domain
    #      ... o xxxxxxxxxxxxx|xxxxxxxxxxxxx||        ||xxxxxxxxxxxxx|xxxxxxxxxxxxx o ...
    #                  |      |      A                        A      |      |
    #                  |      |      |                        |      |      |
    #              receiving  |   sending                  sending   |  receiving
    #                  |      |      |                        |      |      |
    #                  V      |      |                        |      |      V
    #          ||============================= PROC I =============================||
    #          ||<----halo--->|<---------------domain--------------->|<---halo---->||
    #          ||xxxxxxxxxxxxx|xxxxxxxxxxxxx o ...... o xxxxxxxxxxxxx|xxxxxxxxxxxxx||
    #            V             V             V          V             V              V
    # py index   0 ...         hlo ...       2*hlo ...  n ...         n+hlo ...      n+2*hlo
    #
    def swap(self, f, nh, tree):
        """
        Convenience function to swap all directions including the corners in a single call.
        A call to this function is equal to calling swapXc, swapYc and swapZc (in this order)
        separately.

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """
        # 
        self.swapXc(f, nh, tree)
        ndim = tree['eqns']['ndim']
        if ndim > 1:
            self.swapYc(f, nh, tree)
        if ndim > 2:
            self.swapZc(f, nh, tree)

    def swapX(self, f, nh, tree):
        """
        Swaps array "f(X,:)" (a function of position vector X ) in the x-direction
        over its first ":" variables with halo size "nh".

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """
        wp = tree['misc']['working precision']
        ndim = tree['eqns']['ndim']
        iSwap = tree['mpi']['dMpi'].iSwap

        nx = self.nx
        ny = self.ny
        nz = self.nz

        if ndim == 1:
            nv = np.size(f, 1)
            # positive x-dir
            buf = np.empty(nh, dtype=wp)
            buf = f[nx:nh + nx, 0:nv].copy().reshape((nh * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxp,
                    sendtag=0,
                    source=self.neighxm,
                    recvtag=0,
                    status=None)
            if iSwap['i1']:
                f[0:nh, 0:nv] = buf.reshape((nh, nv)).copy()
            # negative x-dir
            buf = np.empty(nh, dtype=wp)
            buf = f[nh:2 * nh, 0:nv].copy().reshape((nh * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxm,
                    sendtag=0,
                    source=self.neighxp,
                    recvtag=0,
                    status=None)
            if iSwap['imax']:
                f[nx + nh:nx + 2 * nh, 0:nv] = buf.reshape((nh, nv)).copy()
        elif ndim == 2:
            nv = np.size(f, 2)
            # positive x-dir
            buf = np.empty(nh * ny * nv, dtype=wp)
            buf = f[nx:nh + nx, nh:nh + ny,
                    0:nv].copy().reshape((nh * ny * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxp,
                    sendtag=0,
                    source=self.neighxm,
                    recvtag=0,
                    status=None)
            if iSwap['i1']:
                f[0:nh, nh:nh + ny, 0:nv] = buf.reshape((nh, ny, nv)).copy()
            # negative x-dir
            buf = np.empty(nh * ny * nv, dtype=wp)
            buf = f[nh:2 * nh, nh:nh + ny,
                    0:nv].copy().reshape((nh * ny * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxm,
                    sendtag=0,
                    source=self.neighxp,
                    recvtag=0,
                    status=None)
            if iSwap['imax']:
                f[nx + nh:nx + 2 * nh, nh:nh + ny,
                    0:nv] = buf.reshape((nh, ny, nv)).copy()
        elif ndim == 3:
            nv = np.size(f, 3)
            # positive x-dir
            buf = np.empty(nh * ny * nz * nv, dtype=wp)
            # t0=timer()
            # buf = f[nx:nh+nx,nh:nh+ny,nh:nh+nz,0:nv].copy().reshape((nh*ny*nz*nv,1))
            # print("f is ",f[:,:,:,0])
            # sys.stdout.flush()
            # sys.exit()
            dnamiF.pack(
                buf,
                f,
                nx,
                nh + nx,
                nh,
                nh + ny,
                nh,
                nh + nz,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
            # t1=timer()
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxp,
                    sendtag=0,
                    source=self.neighxm,
                    recvtag=0,
                    status=None)
            # t2=timer()
            # f[0:nh,nh:nh+ny,nh:nh+nz,0:nv] = buf.reshape((nh,ny,nz,nv)).copy()
            if iSwap['i1']:
                dnamiF.unpack(
                    buf,
                    f,
                    0,
                    nh,
                    nh,
                    nh + ny,
                    nh,
                    nh + nz,
                    nx + 2 * nh,
                    ny + 2 * nh,
                    nz + 2 * nh,
                    nv)
            # print("f is ",f[:,:,:,0])
            # sys.stdout.flush()
            # sys.exit()
            # t3=timer()
            # negative x-dir
            buf = np.empty(nh * ny * nz * nv, dtype=wp)
            # t4=timer()
            # buf = f[nh:2*nh,nh:nh+ny,nh:nh+nz,0:nv].copy().reshape((nh*ny*nz*nv,1))
            dnamiF.pack(
                buf,
                f,
                nh,
                2 * nh,
                nh,
                nh + ny,
                nh,
                nh + nz,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
            # t5=timer()
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxm,
                    sendtag=0,
                    source=self.neighxp,
                    recvtag=0,
                    status=None)
            # t6=timer()
            # f[nx+nh:nx+2*nh,nh:nh+ny,nh:nh+nz,0:nv] = buf.reshape((nh,ny,nz,nv)).copy()
            if iSwap['imax']:
                dnamiF.unpack(
                    buf,
                    f,
                    nx + nh,
                    nx + 2 * nh,
                    nh,
                    nh + ny,
                    nh,
                    nh + nz,
                    nx + 2 * nh,
                    ny + 2 * nh,
                    nz + 2 * nh,
                    nv)
            # t7=timer()
            # print("in mpi x",(t1-t0+t3-t2+t5-t4+t7-t6)/(128**3))
        else:
            if self.ioproc:
                print('[error] (swapX) ndim (', ndim, ') should be 1, 2 or 3')
            sys.exit()

    def swapY(self, f, nh, tree):
        """
        Swaps the array "f(X,:)" (a function of position vector X ) in the y-direction
        over its first ":" variables with halo size "nh". This function does not swap
        the corners.

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """
        #               ->                                   ->
        wp = tree['misc']['working precision']
        ndim = tree['eqns']['ndim']
        iSwap = tree['mpi']['dMpi'].iSwap

        nx = self.nx
        ny = self.ny
        nz = self.nz
        if ndim == 2:
            nv = np.size(f, 2)
            # positive y-dir
            buf = np.empty(nh * nx * nv, dtype=wp)
            buf = f[nh:nh + nx, ny:nh + ny,
                    0:nv].copy().reshape((nh * nx * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighyp,
                    sendtag=0,
                    source=self.neighym,
                    recvtag=0,
                    status=None)
            if iSwap['j1']:
                f[nh:nh + nx, 0:nh, 0:nv] = buf.reshape((nx, nh, nv)).copy()
            # negative y-dir
            buf = np.empty(nh * nx * nv, dtype=wp)
            buf = f[nh:nh + nx, nh:2 * nh,
                    0:nv].copy().reshape((nh * nx * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighym,
                    sendtag=0,
                    source=self.neighyp,
                    recvtag=0,
                    status=None)
            if iSwap['jmax']:
                f[nh:nh + nx, ny + nh:ny + 2 * nh,
                    0:nv] = buf.reshape((nx, nh, nv)).copy()
        elif ndim == 3:
            nv = np.size(f, 3)
            # positive y-dir
            buf = np.empty(nh * nx * nz * nv, dtype=wp)
            buf = f[nh:nh + nx, ny:nh + ny, nh:nh + nz,
                    0:nv].copy().reshape((nh * nx * nz * nv, 1))
            # dnamiF.pack(buf,f,nh,nh+nx,ny,nh+ny,nh,nh+nz,nx+2*nh,ny+2*nh,nz+2*nh,nv)
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighyp,
                    sendtag=0,
                    source=self.neighym,
                    recvtag=0,
                    status=None)
            if iSwap['j1']:
                f[nh:nh + nx, 0:nh, nh:nh + nz,
                    0:nv] = buf.reshape((nx, nh, nz, nv)).copy()
            # dnamiF.unpack(buf,f,nh,nh+nx,0,nh,nh,nh+nz,nx+2*nh,ny+2*nh,nz+2*nh,nv)
            # negative y-dir
            buf = np.empty(nh * nx * nz * nv, dtype=wp)
            buf = f[nh:nh + nx, nh:2 * nh, nh:nh + nz,
                    0:nv].copy().reshape((nh * nx * nz * nv, 1))
            # dnamiF.pack(buf,f,nh,nh+nx,nh,2*nh,nh,nh+nz,nx+2*nh,ny+2*nh,nz+2*nh,nv)
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighym,
                    sendtag=0,
                    source=self.neighyp,
                    recvtag=0,
                    status=None)
            if iSwap['jmax']:
                f[nh:nh + nx, ny + nh:ny + 2 * nh, nh:nh + nz,
                    0:nv] = buf.reshape((nx, nh, nz, nv)).copy()
                # dnamiF.unpack(buf,f,nh,nh+nx,ny+nh,ny+2*nh,nh,nh+nz,nx+2*nh,ny+2*nh,nz+2*nh,nv)

    def swapZ(self, f, nh, tree):
        """
        Swaps array "f(X,:)" (a function of position vector X ) in the z-direction
        over its first ":" variables with halo size "nh".
        This function does not swap the corners.

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """
        wp = tree['misc']['working precision']
        nx = self.nx
        ny = self.ny
        nz = self.nz
        iSwap = tree['mpi']['dMpi'].iSwap

        nv = np.size(f, 3)
        # positive z-dir
        buf = np.empty(nh * nx * ny * nv, dtype=wp)
        buf = f[nh:nh + nx, nh:nh + ny, nz:nh + nz,
                0:nv].copy().reshape((nh * nx * ny * nv, 1))
        if self.iMpi:
            self.comm_torus.Sendrecv_replace(
                buf,
                self.neighzp,
                sendtag=0,
                source=self.neighzm,
                recvtag=0,
                status=None)
        if iSwap['k1']:
            f[nh:nh + nx, nh:nh + ny, 0:nh,
                0:nv] = buf.reshape((nx, ny, nh, nv)).copy()
        # negative z-dir
        buf = np.empty(nh * nx * ny * nv, dtype=wp)
        buf = f[nh:nh + nx, nh:nh + ny, nh:2 * nh,
                0:nv].copy().reshape((nh * nx * ny * nv, 1))
        if self.iMpi:
            self.comm_torus.Sendrecv_replace(
                buf,
                self.neighzm,
                sendtag=0,
                source=self.neighzp,
                recvtag=0,
                status=None)
        if iSwap['kmax']:
            f[nh:nh + nx, nh:nh + ny, nz + nh:nz + 2 * nh,
                0:nv] = buf.reshape((nx, ny, nh, nv)).copy()

# Swap with corners:

    def swapXc(self, f, nh, tree):
        """
        Swaps array "f(X,:)" (a function of position vector X ) in the x-direction
        over its first ":" variables with halo size "nh". This swap function will
        swap the corners.

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """

        wp = tree['misc']['working precision']
        ndim = tree['eqns']['ndim']
        iSwap = tree['mpi']['dMpi'].iSwap

        nx = self.nx
        ny = self.ny
        nz = self.nz

        if ndim == 1:
            nv = np.size(f, 1)
            # positive x-dir
            buf = np.empty(nh, dtype=wp)
            buf = f[nx:nh + nx, 0:nv].copy().reshape((nh * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxp,
                    sendtag=0,
                    source=self.neighxm,
                    recvtag=0,
                    status=None)
            if iSwap['i1']:
                f[0:nh, 0:nv] = buf.reshape((nh, nv)).copy()
            # negative x-dir
            buf = np.empty(nh, dtype=wp)
            buf = f[nh:2 * nh, 0:nv].copy().reshape((nh * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxm,
                    sendtag=0,
                    source=self.neighxp,
                    recvtag=0,
                    status=None)
            if iSwap['imax']:
                f[nx + nh:nx + 2 * nh, 0:nv] = buf.reshape((nh, nv)).copy()
        elif ndim == 2:
            nv = np.size(f, 2)
            # positive x-dir
            buf = np.empty(nh * (ny + 2 * nh) * nv, dtype=wp)
            buf = f[nx:nh + nx, 0:2 * nh + ny,
                    0:nv].copy().reshape((nh * (ny + 2 * nh) * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxp,
                    sendtag=0,
                    source=self.neighxm,
                    recvtag=0,
                    status=None)
            if iSwap['i1']:
                f[0:nh, 0:2 * nh + ny,
                    0:nv] = buf.reshape((nh, ny + 2 * nh, nv)).copy()
            # negative x-dir
            buf = np.empty(nh * (ny + 2 * nh) * nv, dtype=wp)
            buf = f[nh:2 * nh, 0:2 * nh + ny,
                    0:nv].copy().reshape((nh * (ny + 2 * nh) * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxm,
                    sendtag=0,
                    source=self.neighxp,
                    recvtag=0,
                    status=None)
            if iSwap['imax']:
                f[nx + nh:nx + 2 * nh, 0:2 * nh + ny,
                    0:nv] = buf.reshape((nh, ny + 2 * nh, nv)).copy()
        elif ndim == 3:
            nv = np.size(f, 3)
            # positive x-dir
            buf = np.empty(nh * (ny + 2 * nh) * nz * nv, dtype=wp)
            # t0=timer()
            # buf = f[nx:nh+nx,0:2*nh+ny,nh:nh+nz,0:nv].copy().reshape((nh*(ny+2*nh)*nz*nv,1))
            # print("f is ",f[:,:,:,0])
            # sys.stdout.flush()
            # sys.exit()
            dnamiF.pack(
                buf,
                f,
                nx,
                nh + nx,
                0,
                2 * nh + ny,
                nh,
                nh + nz,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
            # t1=timer()
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxp,
                    sendtag=0,
                    source=self.neighxm,
                    recvtag=0,
                    status=None)
            # t2=timer()
            if iSwap['i1']:
                # f[0:nh,0:2*nh+ny,nh:nh+nz,0:nv] = buf.reshape((nh,ny+2*nh,nz,nv)).copy()
                dnamiF.unpack(
                    buf,
                    f,
                    0,
                    nh,
                    0,
                    2 * nh + ny,
                    nh,
                    nh + nz,
                    nx + 2 * nh,
                    ny + 2 * nh,
                    nz + 2 * nh,
                    nv)
            # print("f is ",f[:,:,:,0])
            # sys.stdout.flush()
            # sys.exit()
            # t3=timer()
            # negative x-dir
            buf = np.empty(nh * (ny + 2 * nh) * nz * nv, dtype=wp)
            # t4=timer()
            # buf = f[nh:2*nh,0:2*nh+ny,nh:nh+nz,0:nv].copy().reshape((nh*(ny+2*nh)*nz*nv,1))
            dnamiF.pack(
                buf,
                f,
                nh,
                2 * nh,
                0,
                2 * nh + ny,
                nh,
                nh + nz,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
            # t5=timer()
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighxm,
                    sendtag=0,
                    source=self.neighxp,
                    recvtag=0,
                    status=None)
            # t6=timer()
            if iSwap['imax']:
                # f[nx+nh:nx+2*nh,0:2*nh+ny,nh:nh+nz,0:nv] = buf.reshape((nh,ny+2*nh,nz,nv)).copy()
                dnamiF.unpack(
                    buf,
                    f,
                    nx + nh,
                    nx + 2 * nh,
                    0,
                    2 * nh + ny,
                    nh,
                    nh + nz,
                    nx + 2 * nh,
                    ny + 2 * nh,
                    nz + 2 * nh,
                    nv)
            # t7=timer()
            # print("in mpi x",(t1-t0+t3-t2+t5-t4+t7-t6)/(128**3))
        else:
            if self.ioproc:
                print('[error] (swapX) ndim (', ndim, ') should be 1, 2 or 3')
            sys.exit()

    def swapYc(self, f, nh, tree):
        """
        Swaps array "f(X,:)" (a function of position vector X ) in the y-direction
        over its first ":" variables with halo size "nh". This swap function will
        swap the corners.

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """
        wp = tree['misc']['working precision']
        ndim = tree['eqns']['ndim']
        iSwap = tree['mpi']['dMpi'].iSwap

        nx = self.nx
        ny = self.ny
        nz = self.nz
        if ndim == 2:
            nv = np.size(f, 2)
            # positive y-dir
            buf = np.empty(nh * (nx + 2 * nh) * nv, dtype=wp)
            buf = f[0:2 * nh + nx, ny:nh + ny,
                    0:nv].copy().reshape((nh * (nx + 2 * nh) * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighyp,
                    sendtag=0,
                    source=self.neighym,
                    recvtag=0,
                    status=None)
            if iSwap['j1']:
                f[0:2 * nh + nx, 0:nh,
                    0:nv] = buf.reshape((nx + 2 * nh, nh, nv)).copy()
            # negative y-dir
            buf = np.empty(nh * (nx + 2 * nh) * nv, dtype=wp)
            buf = f[0:2 * nh + nx, nh:2 * nh,
                    0:nv].copy().reshape((nh * (nx + 2 * nh) * nv, 1))
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighym,
                    sendtag=0,
                    source=self.neighyp,
                    recvtag=0,
                    status=None)
            if iSwap['jmax']:
                f[0:2 * nh + nx, ny + nh:ny + 2 * nh,
                    0:nv] = buf.reshape((nx + 2 * nh, nh, nv)).copy()
        elif ndim == 3:
            nv = np.size(f, 3)
            # positive y-dir
            buf = np.empty(nh * (nx + 2 * nh) * nz * nv, dtype=wp)
            # t0=timer()
            # buf = f[0:2*nh+nx,ny:nh+ny,nh:nh+nz,0:nv].copy().reshape((nh*(nx+2*nh)*nz*nv,1))
            dnamiF.pack(
                buf,
                f,
                0,
                2 * nh + nx,
                ny,
                nh + ny,
                nh,
                nh + nz,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
            # t1=timer()
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighyp,
                    sendtag=0,
                    source=self.neighym,
                    recvtag=0,
                    status=None)
            # t2=timer()
            if iSwap['j1']:
                # f[0:2*nh+nx,0:nh,nh:nh+nz,0:nv] = buf.reshape((nx+2*nh,nh,nz,nv)).copy()
                dnamiF.unpack(
                    buf,
                    f,
                    0,
                    2 * nh + nx,
                    0,
                    nh,
                    nh,
                    nh + nz,
                    nx + 2 * nh,
                    ny + 2 * nh,
                    nz + 2 * nh,
                    nv)
            # t3=timer()
            # negative y-dir
            buf = np.empty(nh * (nx + 2 * nh) * nz * nv, dtype=wp)
            # t4=timer()
            # buf = f[0:2*nh+nx,nh:2*nh,nh:nh+nz,0:nv].copy().reshape((nh*(nx+2*nh)*nz*nv,1))
            dnamiF.pack(
                buf,
                f,
                0,
                2 * nh + nx,
                nh,
                2 * nh,
                nh,
                nh + nz,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
            # t5=timer()
            if self.iMpi:
                self.comm_torus.Sendrecv_replace(
                    buf,
                    self.neighym,
                    sendtag=0,
                    source=self.neighyp,
                    recvtag=0,
                    status=None)
            # t6=timer()
            if iSwap['jmax']:
                dnamiF.unpack(
                    buf,
                    f,
                    0,
                    2 * nh + nx,
                    ny + nh,
                    ny + 2 * nh,
                    nh,
                    nh + nz,
                    nx + 2 * nh,
                    ny + 2 * nh,
                    nz + 2 * nh,
                    nv)
                # f[0:2*nh+nx,ny+nh:ny+2*nh,nh:nh+nz,0:nv] = buf.reshape((nx+2*nh,nh,nz,nv)).copy()
            # t7=timer()
            # print("in mpi y",(t1-t0+t3-t2+t5-t4+t7-t6)/(128**3))

    def swapZc(self, f, nh, tree):
        """
        Swaps array "f(X,:)" (a function of position vector X ) in the z-direction
        over its first ":" variables with halo size "nh". This swap function will
        swap the corners.

        Args:
          f: The input and output vector for the halo data swap
          nh: The width of the halo data regions
          tree: The dnami tree data structure
        """
        wp = tree['misc']['working precision']
        iSwap = tree['mpi']['dMpi'].iSwap

        nx = self.nx
        ny = self.ny
        nz = self.nz
        nv = np.size(f, 3)
        # positive z-dir
        buf = np.empty(nh * (nx + 2 * nh) * (ny + 2 * nh) * nv, dtype=wp)
        # t0=timer()
        # buf = f[0:2*nh+nx,0:2*nh+ny,nz:nh+nz,0:nv].copy().reshape((nh*(nx+2*nh)*(ny+2*nh)*nv,1))
        dnamiF.pack(
            buf,
            f,
            0,
            2 * nh + nx,
            0,
            2 * nh + ny,
            nz,
            nh + nz,
            nx + 2 * nh,
            ny + 2 * nh,
            nz + 2 * nh,
            nv)
        # t1=timer()
        if self.iMpi:
            self.comm_torus.Sendrecv_replace(
                buf,
                self.neighzp,
                sendtag=0,
                source=self.neighzm,
                recvtag=0,
                status=None)
        # t2=timer()
        if iSwap['k1']:
            # f[0:2*nh+nx,0:2*nh+ny,0:nh,0:nv] = buf.reshape((nx+2*nh,ny+2*nh,nh,nv)).copy()
            dnamiF.unpack(
                buf,
                f,
                0,
                2 * nh + nx,
                0,
                2 * nh + ny,
                0,
                nh,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)
        # t3=timer()
        # negative z-dir
        buf = np.empty(nh * (nx + 2 * nh) * (ny + 2 * nh) * nv, dtype=wp)
        # t4=timer()
        # buf = f[0:2*nh+nx,0:2*nh+ny,nh:2*nh,0:nv].copy().reshape((nh*(nx+2*nh)*(ny+2*nh)*nv,1))
        dnamiF.pack(
            buf,
            f,
            0,
            2 * nh + nx,
            0,
            2 * nh + ny,
            nh,
            2 * nh,
            nx + 2 * nh,
            ny + 2 * nh,
            nz + 2 * nh,
            nv)
        # t5=timer()
        if self.iMpi:
            self.comm_torus.Sendrecv_replace(
                buf,
                self.neighzm,
                sendtag=0,
                source=self.neighzp,
                recvtag=0,
                status=None)
        # t6=timer()
        if iSwap['kmax']:
            # f[0:2*nh+nx,0:2*nh+ny,nz+nh:nz+2*nh,0:nv] = buf.reshape((nx+2*nh,ny+2*nh,nh,nv)).copy()
            dnamiF.unpack(
                buf,
                f,
                0,
                2 * nh + nx,
                0,
                2 * nh + ny,
                nz + nh,
                nz + 2 * nh,
                nx + 2 * nh,
                ny + 2 * nh,
                nz + 2 * nh,
                nv)


            # t7=timer()
            # print("in mpi z",(t1-t0+t3-t2+t5-t4+t7-t6)/(128**3))
