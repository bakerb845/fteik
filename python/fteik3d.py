#!/usr/bin/env python3
"""
Python interface to the 3D fteik eikonal equation Fortran library.

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
                                                     c_double,  #convTol
                                                     c_int,     #verbosity
                                                     POINTER(c_int) #ierr
                                                    )
        lib.fteik_solver3d_setReceivers64f.argtypes = (c_int,             #nrec
                                                       POINTER(c_double), #zrec
                                                       POINTER(c_double), #xrec
                                                       POINTER(c_double), #yrec
                                                       POINTER(c_int)     #ierr
                                                      )
        lib.fteik_solver3d_solveSourceLSM.argtypes = (c_int,         #(Fortran) source number
                                                      POINTER(c_int) #ierr
                                                     )
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
        self.nrec = 0
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
        self.nrec = 0
        self.linit = False
        self.fteik3d.fteik_solver3d_free()
        return

    def initialize(self, nx, ny, nz, dx, dy, dz,
                   x0 = 0.0, y0 = 0.0, z0 = 0.0,
                   nsweep = 2, eps = 5.0, convTol = 0.0, 
                   verbose=0):
        """
        Initializes the 3D eikonal solver.  This will trigger a graph ordering 
        operation which can be fairly expensive and should only be done once.

        Required Arguments 
        nx : int
           Number of x grid points in the travel time field.
           This must be at least 3.
        ny : int
           Number of y grid points in the travel time field.
           This must be at least 3.
        nz : int
           Number of z grid points in the travel time field. 
           This must be at least 3. 
        dx : float
           Grid spacing (meters) in x.  This cannot be 0.0.
        dy : float
           Grid spacing (meters) in y.  This cannot be 0.0
        dz : float
           Grid spacing (meters) in z.  This cannot be 0.0.

        Optional Arguments 
        x0 : float
           Model origin in x (meters). 
        y0 : float
           Model origin in y (meters). 
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
        convTol : float
           The Gauss-Seidel method will terminate if the updates are less than
           convTol (seconds).
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
                                                  convTol, verbose,
                                                  byref(ierr))
        if (ierr.value != 0):
            print("Error initializing solver")
            return -1
        self.nx = nx
        self.ny = ny
        self.nz = nz
        return 0

    def setVelocityModel(self, vel):
        """
        Sets the nodal or cell-based velocity model.  
        Input
        -----
        vel : numpy matrix
           This is the [nz-1 x ny-1 x nx-1] cell-based velocity model or
           [nz x ny x nx] node-based velocity model in meters/second.
           Note, that if a nodal model is specified then the solver
           will perform a simple averaging scheme to convert it to a
           cell-based model. 
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

        Returns
        -------
        errorAll : int
           Error counter where 0 indicates a success.
        """ 
        ierr = c_int(1)
        errorAll = 0
        for i in range(self.nsrc):
            isrc = i + 1
            self.fteik3d.fteik_solver3d_solveSourceLSM(isrc, ierr)
            if (ierr.value != 0):
                print("Failed to solve for source %d"%i+1)
                errorAll = errorAll + 1
        return errorAll

    def getTravelTimeField(self):
        """
        Extracts the travel time field from the solver.

        Result
        ttimes : numpy matrix
            A [nz x ny x nx] matrix containing the travel times (seconds)
            from the source to each node in the model.
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

    def getTravelTimes(self):
        """
        Gets the travel times at the receivers.

        Results
        ttr : array_like
           On successful exit this is the travel times (seconds) from the source
           to the receiver.
        """ 
        nrec = self.nrec
        if (nrec < 1):
            print("No receiver locations set")
            return None
        ttr = ascontiguousarray(zeros(nrec), dtype='float64')
        ttrPointer = ttr.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.fteik3d.fteik_solver3d_getTravelTimes64f(nrec, ttrPointer, ierr)
        if (ierr.value != 0):
            print("Error getting travel times")
            return None
        return ttr

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
        nsrc = min(len(xsrc), len(ysrc), len(zsrc))
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
            return -1
        self.nsrc = nsrc
        return 0

    def setReceivers(self, xrec, yrec, zrec):
        """
        Sets the receiver positions.

        Input
        xrec : array_like
          x receiver locations (meters).
        yrec : array_like
          y receiver locations (meters).
        zrec : array_like
          z receiver locations (meters).

        Returns
        ierr : int 
           0 indicates success.
        """
        xrec = ascontiguousarray(xrec, float64)
        yrec = ascontiguousarray(yrec, float64)
        zrec = ascontiguousarray(zrec, float64)
        nrec = min(len(xrec), len(yrec), len(zrec))
        if (len(xrec) != len(zrec)):
            print("Inconsistent array lengths")
        xrecPointer = xrec.ctypes.data_as(POINTER(c_double))
        yrecPointer = yrec.ctypes.data_as(POINTER(c_double))
        zrecPointer = zrec.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.fteik3d.fteik_solver3d_setReceivers64f(nrec,
                                                    zrecPointer,
                                                    xrecPointer,
                                                    yrecPointer,
                                                    ierr)
        if (ierr.value != 0):
            print("Error setting receivers")
            return -1
        self.nrec = nrec
        return 0

