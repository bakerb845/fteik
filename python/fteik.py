#!/usr/bin/python3
"""
Python interface to the fteik eikonal equation Fortran library.

Copyright: Ben Baker distributed under the MIT license.
"""
import ctypes
from ctypes import c_int
from ctypes import c_double
from ctypes import byref
from ctypes import POINTER
from numpy import zeros
from numpy import array
from numpy import reshape
from numpy import float64
from numpy import ascontiguousarray 
from numpy import linspace
import os

class fteik2d:
    def __init__(self,
                 fteik_library='/usr/local/lib/libfteik_shared.so'):
        lib = ctypes.cdll.LoadLibrary(fteik_library)
        # Make the interfaces
        lib.fteik_solver2d_initialize64f.argtypes = (c_int,     #nz
                                                     c_int,     #nx
                                                     c_double,  #z0
                                                     c_double,  #x0
                                                     c_double,  #dz
                                                     c_double,  #dx
                                                     c_int,     #nsweeps
                                                     c_double,  #eps
                                                     c_int,     #verbosity
                                                     POINTER(c_int) #ierr
                                                    )
        #lib.fteik_solver2d_setVelocityModel64f.argtypes = (c_int,             #ncell
        #                                                   POINTER(c_double), #velocity
        #                                                   POINTER(c_int)     #ierr
        #                                                  )
        lib.fteik_solver2d_setCellVelocityModel64f.argtypes = (c_int,             #ncell
                                                               c_int,             #order
                                                               POINTER(c_double), #velocity
                                                               POINTER(c_int)     #ierr
                                                              )
        lib.fteik_solver2d_setNodalVelocityModel64f.argtypes = (c_int,             #ngrd
                                                                c_int,             #order
                                                                POINTER(c_double), #velocity
                                                                POINTER(c_int)     #ierr
                                                               )
        lib.fteik_solver2d_setReceivers64f.argtypes = (c_int,             #nrec
                                                       POINTER(c_double), #zrec
                                                       POINTER(c_double), #xrec
                                                       POINTER(c_int)     #ierr
                                                      )
        lib.fteik_solver2d_solveSourceLSM.argtypes = (c_int,         #(Fortran) source number
                                                      POINTER(c_int) #ierr
                                                     )
        lib.fteik_solver2d_solveSourceFSM.argtypes = (c_int,         #(Fortran) source number
                                                      POINTER(c_int) #ierr
                                                     )
        lib.fteik_solver2d_getTravelTimeField64f.argtypes = (c_int,             #ngrd
                                                             c_int,             #order
                                                             POINTER(c_double), #ttimes
                                                             POINTER(c_int)     #ierr
                                                            )
        lib.fteik_solver2d_setSources64f.argtypes = (c_int,             #nsrc
                                                     POINTER(c_double), #zsrc
                                                     POINTER(c_double), #xsrc
                                                     POINTER(c_int)     #ierr
                                                    )
        lib.fteik_solver2d_getTravelTimes64f.argtypes = (c_int,             #nrec
                                                         POINTER(c_double), #ttr
                                                         POINTER(c_int)     #ierr
                                                        )
        lib.fteik_solver2d_getNumberOfReceivers.argtypes = (POINTER(c_int), #nrec 
                                                            POINTER(c_int)  #ierr
                                                           )
        lib.fteik_solver2d_free.argtypes = None
        self.linit = False
        self.nx = 0
        self.nz = 0
        self.nsrc = 0
        self.fteik2d = lib
        return

    def __enter__(self):
        return self

    def __exit__(self):
        self.fteik2d.fteik_solver2d_free() 
        return
    #####################################################################################
    #                                The Python Interfaces                              #
    #####################################################################################
    def free(self):
        """
        Frees the internal Fortran module.
        """
        self.linit = False
        self.nx = 0
        self.nz = 0
        self.nsrc = 0
        self.fteik2d.fteik_solver2d_free()
        return

    def initialize(self, nx, nz, dx, dz,
                   x0 = 0.0, z0 = 0.0,
                   nsweep = 2, eps = 3.0, verbose=0):
        """
        Initializes the 2D solver eikonal solver.

        Required Arguments 
        nx : int
           Number of x grid points in the travel time field.
           This must be at least 3.
        nz : int
           Number of z grid points in the travel time field. 
           This must be at least 3. 
        dx : float
           Grid spacing (meters) in x.  This cannot be 0.0.
        dz : float
           Grid spacing (meters) in z.  This cannot be 0.0.

        Optional Arguments 
        x0 : float
           Model origin in x (meters). 
        z0 : float
           Model origin in z (meters).
        nsweep : int
           Number of Gauss-Seidel iterations.
        eps : float
           Number of grid-points around source where, on initialization, the
           the solver switches from a spherical solver to a cartesian solver.
           This variable can be thought of in the sense that the eikonal 
           equation assumes a plane wave solution to the acoustic wave equation.
           In the source region there will be high-curvature which is difficult
           to approximate in cartesian coordinates.

        verbose : int
           Controls verbosity where 0 is quiet and 1, 2, ... correspond to
           an increasing number of messages from the module to standard out.
        """
        ierr = c_int(1)
        self.fteik2d.fteik_solver2d_initialize64f(nz, nx,
                                                  z0, x0,
                                                  dz, dx, 
                                                  nsweep, eps,
                                                  verbose, byref(ierr))
        if (ierr.value != 0):
            print("Error initializing solver")
            return -1 
        self.nx = nx
        self.nz = nz
        self.nsrc = 0
        return 0

    def setVelocityModel(self, vel):
        """
        Sets the [nz-1 x nx-1] cell-based model velocity model or the [nz x nx]
        nodal-based model.

        Input
        -----
        vel : np.matrix
           [nz-1 x nx-1] cell-based model or [nz x nx] nodal-based model of
           velocities (m/s).

        Returns
        -------
        ierr : int
           0 indicates success.
        """
        vel = reshape(vel, vel.size, order='F')   # Make it a 1D [nz-1 x nx-1] array
        vel = ascontiguousarray(vel, float64)
        lhaveGrid = False
        ngrd = len(vel)
        ncell = len(vel)
        if (ncell != (self.nx - 1)*(self.nz - 1)):
            if (ngrd  != self.nx*self.nz):
                print("Expecing %d grid points or %d cells"%
                      self.nx*self.ny, (self.nx - 1)*(self.nz - 1))
            else:
                lhaveGrid = True
        velPointer = vel.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        order = order = int(0) # FTEIK_ZX_ORDERING
        if (lhaveGrid):
            self.fteik2d.fteik_solver2d_setNodalVelocityModel64f(ngrd, order,
                                                                 velPointer,
                                                                 byref(ierr)) 
        else:
            self.fteik2d.fteik_solver2d_setCellVelocityModel64f(ncell, order,
                                                                velPointer,
                                                                byref(ierr))
        #self.fteik2d.fteik_solver2d_setVelocityModel64f(ncell, velPointer,
        #                                                byref(ierr))
        if (ierr.value != 0):
            print("Error setting velocity")
            return -1
        return 0

    def solveLSM(self):
        """
        Solves the eikonal equation with the level-set method
        """
        ierr = c_int(1) 
        for i in range(self.nsrc):
            isrc = i + 1
            self.fteik2d.fteik_solver2d_solveSourceLSM(isrc, ierr)

    def getTravelTimeField(self):
        """
        Extracts the travel time field from the solver.
        """
        ngrd = self.nx*self.nz
        ttimes = ascontiguousarray(zeros(ngrd), dtype='float64')
        ttimesPointer = ttimes.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        order = int(0) # FTEIK_ZX_ORDERING 
        self.fteik2d.fteik_solver2d_getTravelTimeField64f(ngrd, order,
                                                          ttimesPointer, ierr)
        if (ierr.value != 0):
            print("Error getting travel time field")
            return None
        ttimes = reshape(ttimes, [self.nz, self.nx], order='F')
        return ttimes

    def setSources(self, xsrc, zsrc):
        """
        Sets the source locations on model.

        Input
        xsrc : array_like
          x source locations (meters).
        zsrc : array_like
          z source locations (meters).

        Returns
        ierr : int
          0 indicates success. 
        """
        xsrc = ascontiguousarray(xsrc, float64)
        zsrc = ascontiguousarray(zsrc, float64)
        nsrc = len(xsrc)
        if (len(xsrc) != len(zsrc)):
            print("Inconsistent array lengths")
        xsrcPointer = xsrc.ctypes.data_as(POINTER(c_double))
        zsrcPointer = zsrc.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.fteik2d.fteik_solver2d_setSources64f(nsrc,
                                                  zsrcPointer, xsrcPointer,
                                                  ierr)
        if (ierr.value != 0):
            print("Error setting sources")
        self.nsrc = nsrc
        return 0 

    def setReceivers(self, xrec, zrec):
        return 

################################################################################
#                                 3D Class                                     #
################################################################################

class fteik3d:
    def __init__(self,
                 fteik_library='/usr/local/lib/libfteik_shared.so'):
        lib = ctypes.cdll.LoadLibrary(fteik_library)
        # Make the interfaces
        lib.fteik_solver3d_initialize64f.argtypes = (c_int,     #nz 
                                                     c_int,     #nx 
                                                     c_int,     #ny
                                                     c_double,  #z0 
                                                     c_double,  #x0 
                                                     c_double,  #y0
                                                     c_double,  #dz 
                                                     c_double,  #dx 
                                                     c_double,  #dz
                                                     c_int,     #nsweeps
                                                     c_double,  #eps
                                                     c_int,     #verbosity
                                                     POINTER(c_int) #ierr
                                                    )
        #lib.fteik_solver3d_setVelocityModel64f.argtypes = (c_int,             #ncell
        #                                                   POINTER(c_double), #velocity
        #                                                   POINTER(c_int)     #ierr
        #                                                  )
        lib.fteik_solver3d_setReceivers64f.argtypes = (c_int,             #nrec
                                                       POINTER(c_double), #zrec
                                                       POINTER(c_double), #xrec
                                                       POINTER(c_double), #yrec
                                                       POINTER(c_int)     #ierr
                                                      )
        lib.fteik_solver3d_solveSourceLSM.argtypes = (c_int,         #(Fortran) source number
                                                      POINTER(c_int) #ierr
                                                     )
        #lib.fteik_solver2d_solveSourceFSM.argtypes = (c_int,         #(Fortran) source number
        #                                              POINTER(c_int) #ierr
        #                                             )
        lib.fteik_solver3d_getTravelTimeField64f.argtypes = (c_int,             #ngrd
                                                             c_int,             #order
                                                             POINTER(c_double), #ttimes
                                                             POINTER(c_int)     #ierr
                                                            )
        lib.fteik_solver3d_setCellVelocityModel64f.argtypes = (c_int,             #ncell
                                                               c_int,             #order
                                                               POINTER(c_double), #velocity
                                                               POINTER(c_int)     #ierr
                                                              )
        lib.fteik_solver3d_setNodalVelocityModel64f.argtypes = (c_int,             #ngrd
                                                                c_int,             #order
                                                                POINTER(c_double), #velocity
                                                                POINTER(c_int)     #ierr
                                                               )
        lib.fteik_solver3d_setSources64f.argtypes = (c_int,             #nsrc
                                                     POINTER(c_double), #zsrc
                                                     POINTER(c_double), #xsrc
                                                     POINTER(c_double), #ysrc
                                                     POINTER(c_int)     #ierr
                                                    )
        lib.fteik_solver3d_getTravelTimes64f.argtypes = (c_int,             #nrec
                                                         POINTER(c_double), #ttr
                                                         POINTER(c_int)     #ierr
                                                        )
        lib.fteik_solver3d_getNumberOfReceivers.argtypes = (POINTER(c_int), #nrec 
                                                            POINTER(c_int)  #ierr
                                                           )

        lib.fteik_solver3d_free.argtypes = None
        self.linit = False
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.nsrc = 0
        self.fteik3d = lib
        return

    def __enter__(self):
        return self

    def __exit__(self):
        self.fteik3d.fteik_solver3d_free()
        return

    def free(self):
        """
        Frees the internal Fortran module.
        """
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.nsrc = 0
        self.linit = False
        self.fteik3d.fteik_solver3d_free()
        return

    def initialize(self, nx, ny, nz, dx, dy, dz,
                   x0 = 0.0, y0 = 0.0, z0 = 0.0,
                   nsweep = 2, eps = 3.0, verbose=0):
        """
        Initializes the 3D eikonal solver.  This will trigger a graph ordering 
        operation which can be fairly expensive and should only be done once.

        Required Arguments 
        nx : int
           Number of x grid points in the travel time field.
           This must be at least 3.
        nz : int
           Number of z grid points in the travel time field. 
           This must be at least 3. 
        dx : float
           Grid spacing (meters) in x.  This cannot be 0.0.
        dz : float
           Grid spacing (meters) in z.  This cannot be 0.0.

        Optional Arguments 
        x0 : float
           Model origin in x (meters). 
        z0 : float
           Model origin in z (meters).
        nsweep : int
           Number of Gauss-Seidel iterations.
        eps : float
           Number of grid-points around source where, on initialization, the
           the solver switches from a spherical solver to a cartesian solver.
           This variable can be thought of in the sense that the eikonal 
           equation assumes a plane wave solution to the acoustic wave equation.
           In the source region there will be high-curvature which is difficult
           to approximate in cartesian coordinates.

        verbose : int
           Controls verbosity where 0 is quiet and 1, 2, ... correspond to
           an increasing number of messages from the module to standard out.

        Returns
        ierr : int
           0 indicates success.
        """
        ierr = c_int(1)
        self.fteik3d.fteik_solver3d_initialize64f(nz, nx, ny,
                                                  z0, x0, y0,
                                                  dz, dx, dy,
                                                  nsweep, eps,
                                                  verbose, byref(ierr))
        if (ierr.value != 0):
            print("Error initializing solver")
            return -1
        self.nx = nx
        self.ny = ny
        self.nz = nz
        return 0

    def setVelocityModel(self, vel):
        """
        Sets the [nz-1 x ny-1 x nx-1] velocity model.
        """
        # Make it a 1D [nz-1 x ny-1 x nx -1] array
        vel = reshape(vel, vel.size, order='F')
        vel = ascontiguousarray(vel, float64)
        lhaveGrid = False
        ngrd = len(vel)
        ncell = len(vel)
        if (ncell != (self.nx - 1)*(self.ny - 1)*(self.nz - 1)):
            if (ngrd  != self.nx*self.ny*self.nz):
                print("Expecing %d grid points or %d cells"%
                      self.nx*self.ny*self.nz,
                      (self.nx - 1)*(self.ny - 1)*(self.nz - 1))
            else:
                lhaveGrid = True
        velPointer = vel.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        order = int(2) # FTEIK_ZYX_ORDERING
        if (lhaveGrid):
            self.fteik3d.fteik_solver3d_setNodalVelocityModel64f(ngrd, order,
                                                                 velPointer,
                                                                 byref(ierr)) 
        else:
            self.fteik3d.fteik_solver3d_setCellVelocityModel64f(ncell, order,
                                                                velPointer,
                                                                byref(ierr))
        if (ierr.value != 0): 
            print("Error setting velocity")
            return -1
        return 0

    def solveLSM(self):
        """
        Solves the eikonal equation with the level-set method
        """
        ierr = c_int(1) 
        for i in range(self.nsrc):
            isrc = i + 1 
            self.fteik3d.fteik_solver3d_solveSourceLSM(isrc, ierr)

    def setSources(self, xsrc, ysrc, zsrc):
        """
        Sets the source locations on model.

        Input
        xsrc : array_like
          x source locations (meters).
        ysrc : array_like
          y source locations (meters).
        zsrc : array_like
          z source locations (meters).

        Returns
        ierr : int
          0 indicates success. 
        """
        xsrc = ascontiguousarray(xsrc, float64)
        ysrc = ascontiguousarray(ysrc, float64)
        zsrc = ascontiguousarray(zsrc, float64)
        nsrc = len(xsrc)
        if (len(xsrc) != len(zsrc)):
            print("Inconsistent array lengths")
        xsrcPointer = xsrc.ctypes.data_as(POINTER(c_double))
        ysrcPointer = ysrc.ctypes.data_as(POINTER(c_double))
        zsrcPointer = zsrc.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.fteik3d.fteik_solver3d_setSources64f(nsrc,
                                                  zsrcPointer,
                                                  xsrcPointer,
                                                  ysrcPointer,
                                                  ierr)
        if (ierr.value != 0):
            print("Error setting sources")
        self.nsrc = nsrc
        return 0

    def getTravelTimeField(self):
        """
        Extracts the travel time field from the solver.
        """
        ngrd = self.nz*self.nx*self.ny
        ttimes = ascontiguousarray(zeros(ngrd), dtype='float64')
        ttimesPointer = ttimes.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        order = int(2) # FTEIK_ZYX_ORDERING 
        self.fteik3d.fteik_solver3d_getTravelTimeField64f(ngrd, order,
                                                          ttimesPointer, ierr)
        if (ierr.value != 0):
            print("Error getting travel time field")
            return None
        ttimes = reshape(ttimes, [self.nz, self.ny, self.nx], order='F')
        return ttimes


if __name__ == "__main__":
    from matplotlib import pyplot as plt
    nx = 101
    nz = 101
    dx = 10.0
    dz = 10.0
    x = linspace(0.0, (nx-1)*dx, nx)
    z = linspace(0.0, (nz-1)*dz, nz)
    xsrc = (nx - 1)*dx/4.0
    zsrc = (nz - 1)*dz/4.0
    vel = zeros([nz,nx]) + 5.e3 #zeros([nz-1,nx-1]) + 5.e3
    fteik2d = fteik2d('/home/bakerb25/C/fteik/lib/libfteik_shared.so')
    fteik2d.initialize(nx, nz, dx, dz, verbose=3)
    fteik2d.setVelocityModel(vel) 
    fteik2d.setSources(xsrc, zsrc)
    fteik2d.solveLSM() 
    ttimes = fteik2d.getTravelTimeField()
    fteik2d.free()
    # plot it
    """
    print(ttimes.shape)
    plt.contourf(x, z, ttimes, 25, cmap=plt.cm.viridis)
    plt.ylim((nz - 1)*dz, 0)
    plt.title('Travel Time Field (s)')
    plt.xlabel('Offset (m)')
    plt.ylabel('Depth (m)')
    plt.colorbar()
    plt.show()
    """

    nx = 51
    ny = 51
    nz = 51
    dx = 10.0
    dy = 10.0
    dz = 10.0
    xsrc = (nx - 1)*dx/4.0
    ysrc = (ny - 1)*dy/2.0
    zsrc = (nz - 1)*dz/4.0
    vmin = 5000.0
    vmax = 6000.0 
    vel = zeros([nz,ny,nx]) + 5.e3 #zeros([nz-1,nx-1,ny-1]) + 5.e3
    for iz in range(nz):
        vel[iz,:,:] = vmin + (vmax - vmin)/((nz-1)*dz)*(iz*dz)
        #print(vmin + (vmax - vmin)/((nz-1)*dz)*(iz*dz))
    fteik3d = fteik3d('/home/bakerb25/C/fteik/lib/libfteik_shared.so')
    fteik3d.initialize(nx, ny, nz, dx, dy, dz, verbose=3) 
    fteik3d.setVelocityModel(vel)
    fteik3d.setSources(xsrc, ysrc, zsrc)
    fteik3d.solveLSM()
    ttimes = fteik3d.getTravelTimeField()
    fteik3d.free()

    """
    print(ttimes.shape)
    xindx = int(xsrc/dx)
    yindx = int(ysrc/dy)
    zindx = int(zsrc/dz)
    x = linspace(0.0, (nx-1)*dx, nx)
    y = linspace(0.0, (ny-1)*dy, ny)
    z = linspace(0.0, (nz-1)*dz, nz)
    # z-x slice
    fig = plt.figure(figsize=[10,10])
    plt.plot([1,1,3])
    plt.subplot(311)
    plt.contourf(x, z, ttimes[:,yindx,:], 25, cmap=plt.cm.viridis)
    plt.ylim((nz - 1)*dz, 0)
    plt.title('X-Z Travel Time Field (s)')
    #plt.xlabel('X-Offset (m)')
    plt.ylabel('Depth (m)')
    plt.colorbar()
    # z-y slice
    plt.subplot(312) 
    plt.contourf(y, z, ttimes[:,:,xindx], 25, cmap=plt.cm.viridis)
    plt.ylim((nz - 1)*dz, 0)
    plt.title('Y-Z Travel Time Field (s)')
    #plt.xlabel('Y-Offset (m)')
    plt.ylabel('Depth (m)')
    plt.colorbar()
    # x-y slice
    plt.subplot(313) 
    plt.contourf(x, y, ttimes[zindx,:,:], 25, cmap=plt.cm.viridis)
    plt.title('X-Y Travel Time Field (s)')
    plt.xlabel('X-Offset (m)')
    plt.ylabel('Y-Offset (m)')
    plt.colorbar()
    plt.show()
    """
