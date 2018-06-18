#!/usr/bin/env python3
"""
Python interface to the 2D fteik eikonal equation Fortran library.

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
                 fteik_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 fteik_library='libfteik_shared.so'):
        # Load the library
        lfound = False
        for path in fteik_path:
            fteik_lib = os.path.join(path, fteik_library)
            if (os.path.isfile(fteik_lib)):
                lfound = True
                break
        if (lfound):
            lib = ctypes.cdll.LoadLibrary(fteik_lib)
        else:
            print("Couldn't find %s"%fteik_library)
            return
        # Make the interfaces
        lib.fteik_solver2d_initialize64f.argtypes = (c_int,     #nz 
                                                     c_int,     #nx 
                                                     c_double,  #z0 
                                                     c_double,  #x0 
                                                     c_double,  #dz 
                                                     c_double,  #dx 
                                                     c_int,     #nsweeps
                                                     c_double,  #eps
                                                     c_double,  #convTol
                                                     c_int,     #verbosity
                                                     POINTER(c_int) #ierr
                                                    )
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
        lib.fteik_solver2d_solveLSM.argtypes = (POINTER(c_int),) # ierr
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
        self.nrec = 0
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
        self.nrec = 0
        self.fteik2d.fteik_solver2d_free()
        return

    def initialize(self, nx, nz, dx, dz,
                   x0 = 0.0, z0 = 0.0,
                   nsweep = 2, eps = 3.0,
                   convTol = 1.e-6, verbose=0):
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
        convTol : float
           The Gauss-Seidel method will terminate if the updates are less than
           convTol (seconds).

        verbose : int
           Controls verbosity where 0 is quiet and 1, 2, ... correspond to
           an increasing number of messages from the module to standard out.
        """
        ierr = c_int(1)
        self.fteik2d.fteik_solver2d_initialize64f(nz, nx,
                                                  z0, x0,
                                                  dz, dx,
                                                  nsweep, eps,
                                                  convTol, verbose, byref(ierr))
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
        self.fteik2d.fteik_solver2d_solveLSM(ierr)
        if (ierr.value != 0):
            print("Error solving eikonal equation")
            return -1
        return 0
        #errorAll = 0
        #for i in range(self.nsrc):
        #    isrc = i + 1
        #    self.fteik2d.fteik_solver2d_solveSourceLSM(isrc, ierr)
        #    if (ierr.value != 0):
        #        print("Failed to solve for source %d"%i+1)
        #        errorAll = errorAll + 1
        #return errorAll

    def getTravelTimeField(self):
        """
        Extracts the travel time field from the solver.

        Result
        ttimes : numpy matrix
           A [nz x nx] matrix containing the travel times (seconds)
           from the source to each node in the model.
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

    def getTravelTimes(self):
        """
        Gets the travel times at the receivers.

        Results
        ttr : matrix.
           On successful exit this is the travel times (seconds) from all
           sources to the receivers.  This is a [nrec x nsrc] matrix.
        """ 
        nrec = self.nrec
        nsrc = self.nsrc
        if (nrec < 1): 
            print("No receiver locations set")
            return None
        if (nsrc < 1):
            print("No sources")
        ttr = ascontiguousarray(zeros(nrec*nsrc), dtype='float64')
        ttrPointer = ttr.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.fteik2d.fteik_solver2d_getTravelTimes64f(nrec, ttrPointer, ierr)
        if (ierr.value != 0): 
            print("Error getting travel times")
            return None
        if (nsrc > 1):
            ttr = reshape(ttr, [self.nrec, self.nsrc], order='F')
        return ttr

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
            return -1
        self.nsrc = nsrc
        return 0

    def setReceivers(self, xrec, zrec):
        """
        Sets the receiver positions.

        Input
        xrec : array_like
           x receiver locations (meters).
        zrec : array_like
           z receiver locations (meters).

        Returns
        ierr : int 
           0 indicates success.
        """
        xrec = ascontiguousarray(xrec, float64)
        zrec = ascontiguousarray(zrec, float64)
        nrec = min(len(xrec), len(zrec))
        if (len(xrec) != len(zrec)):
            print("Inconsistent array lengths")
        xrecPointer = xrec.ctypes.data_as(POINTER(c_double))
        zrecPointer = zrec.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.fteik2d.fteik_solver2d_setReceivers64f(nrec,
                                                    zrecPointer,
                                                    xrecPointer,
                                                    ierr)
        if (ierr.value != 0): 
            print("Error setting receivers")
            return -1
        self.nrec = nrec
        return 0

