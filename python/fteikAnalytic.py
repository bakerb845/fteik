#!/usr/bin/python3
"""
Python interface to the analytic solutions.

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

class fteikAnalytic:
    def __init__(self,
                 fteik_library='/usr/local/lib/libfteik_shared.so'):
        lib = ctypes.cdll.LoadLibrary(fteik_library)
        # Make the interface
        lib.fteik_analytic_initialize64f.argtypes = (c_int,     #nz 
                                                     c_int,     #nx 
                                                     c_int,     #ny 
                                                     c_double,  #z0 
                                                     c_double,  #x0 
                                                     c_double,  #y0 
                                                     c_double,  #dz 
                                                     c_double,  #dx 
                                                     c_double,  #dz 
                                                     c_int,     #verbosity
                                                     POINTER(c_int) #ierr
                                                    ) 
        lib.fteik_analytic_setReceivers64f.argtypes = (c_int,             #nrec
                                                       POINTER(c_double), #zrec
                                                       POINTER(c_double), #xrec
                                                       POINTER(c_double), #yrec
                                                       POINTER(c_int)     #ierr
                                                      )
        lib.fteik_analytic_solveSourceLinearVelocityGradient64f.argtypes = (c_int,    #job
                                                                            POINTER(c_int) #ierr
                                                                           )
        lib.fteik_analytic_solveSourceConstantVelocity64f.argtypes = (c_int,    #job
                                                                      POINTER(c_int) #ierr
                                                                     )
        lib.fteik_analytic_setLinearVelocityGradient64f.argtypes = (c_double,      #v0
                                                                    c_double,      #v1
                                                                    POINTER(c_int) #ierr
                                                                    )
        lib.fteik_analytic_setConstantVelocity64f.argtypes = (c_double,      #v0 
                                                              POINTER(c_int) #ierr
                                                             )
        lib.fteik_analytic_getTravelTimeField64f.argtypes = (c_int,             #ngrd
                                                             c_int,             #order
                                                             POINTER(c_double), #ttimes
                                                             POINTER(c_int)     #ierr
                                                            )
        lib.fteik_analytic_getTravelTimesConstantVel64f.argtypes = (c_int,             #nrec
                                                                    POINTER(c_double), #ttr
                                                                    POINTER(c_int)     #ierr
                                                                   )
        lib.fteik_analytic_getTravelTimesGradientVel64f.argtypes = (c_int,             #nrec
                                                                    POINTER(c_double), #ttr
                                                                    POINTER(c_int)     #ierr
                                                                   )
        lib.fteik_analytic_setSources64f.argtypes = (c_int,             #nsrc
                                                     POINTER(c_double), #zsrc
                                                     POINTER(c_double), #xsrc
                                                     POINTER(c_double), #ysrc
                                                     POINTER(c_int)     #ierr
                                                    )
        lib.fteik_analytic_free.argtypes = None

        self.nsrc = 0
        self.nrec = 0
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.analytic = lib
        return

    def __enter__(self): 
        return self

    def __exit__(self):
        self.analytic.fteik_analytic_free()
        return

    def free(self):
        self.nsrc = 0
        self.nrec = 0
        self.nx = 0
        self.ny = 0
        self.nz = 0
        self.analytic.fteik_analytic_free()
        return

    def initialize(self, nx, ny, nz, dx, dy, dz,
                   x0 = 0.0, y0 = 0.0, z0 = 0.0,
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

        verbose : int
           Controls verbosity where 0 is quiet and 1, 2, ... correspond to
           an increasing number of messages from the module to standard out.

        Returns
        ierr : int
           0 indicates success.
        """
        ierr = c_int(1)
        self.analytic.fteik_analytic_initialize64f(nz, nx, ny,
                                                   z0, x0, y0,
                                                   dz, dx, dy,
                                                   verbose, byref(ierr))
        if (ierr.value != 0):
            print("Error initializing solver")
            return -1
        self.nx = nx
        self.ny = ny
        self.nz = nz
        return 0

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
        self.analytic.fteik_analytic_setSources64f(nsrc,
                                                   zsrcPointer,
                                                   xsrcPointer,
                                                   ysrcPointer,
                                                   ierr)
        if (ierr.value != 0):
            print("Error setting sources")
            return -1
        self.nsrc = nsrc
        return 0

    def setConstantVelocity(self, vel=5000.0):
        """
        Sets the constant velocity model.

        vel : float
           Constant velocity (m/s).
        """
        if (vel <= 0.0):
            print("Velocity must be constant")
        ierr = c_int(1)
        self.analytic.fteik_analytic_setConstantVelocity64f(vel, ierr)
        if (ierr.value != 0): 
            print("Error setting velocity")
            return -1
        return 0

    def setLinearVelocityGradient(self, v1=5000.0, v2=6000.0):
        """ 
        Solves the eikonal equaiton given a linear velocity gradient. 
        """
        if (v1 <= 0.0 or v2 <= 0.0):
            print("Velocities must be positive")
        ierr = c_int(1)
        self.analytic.fteik_analytic_setLinearVelocityGradient64f(v1, v2, ierr)
        if (ierr.value != 0): 
            print("Error setting velocity gradient")
            return -1
        return 0

    def solveConstantVelocity(self):
        """
        Solves the eikonal equation given a contant travel time. 

        Returns
        errorAll : int
           0 indicates success.
        """
        ierr = c_int(1)
        errorAll = 0 
        for i in range(self.nsrc):
            isrc = i + 1 
            self.analytic.fteik_analytic_solveSourceConstantVelocity64f(isrc, ierr)
            if (ierr.value != 0): 
                print("Failed to solve for source %d"%i+1)
                errorAll = errorAll + 1 
        return errorAll

    def solveLinearVelocityGradient(self):
        """
        Solves the eikonal equaiton given a linear velocity gradient. 

        Returns
        errorAll : int
           0 indicates success.
        """
        ierr = c_int(1)
        errorAll = 0
        for i in range(self.nsrc):
            isrc = i + 1 
            self.analytic.fteik_analytic_solveSourceLinearVelocityGradient64f(isrc, ierr)
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
        self.analytic.fteik_analytic_getTravelTimeField64f(ngrd, order,
                                                           ttimesPointer, ierr)
        if (ierr.value != 0): 
            print("Error getting travel time field")
            return None
        if (self.ny > 1):
            ttimes = reshape(ttimes, [self.nz, self.ny, self.nx], order='F')
        else:
            ttimes = reshape(ttimes, [self.nz, self.nx], order='F')
        return ttimes

    def getTravelTimesConstantVelocity(self):
        """
        Gets the travel times at the receivers in constant velocity model.

        Results
        ttr : array_like
           On successful exit this is the travel times (seconds) from the source
           to the receiver.
        """ 
        nsrc = self.nsrc
        nrec = self.nrec
        if (nsrc < 1):
            print("No sources")
            return None
        if (nrec < 1): 
            print("No receiver locations set")
            return None
        ttr = ascontiguousarray(zeros(nrec*nsrc), dtype='float64')
        ttrPointer = ttr.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.analytic.fteik_analytic_getTravelTimesConstantVel64f(nrec, ttrPointer, ierr)
        if (ierr.value != 0): 
            print("Error getting travel times")
            return None
        if (nsrc > 1):
            ttr = reshape(ttr, [self.nrec, self.nsrc], order='F')
        return ttr

    def getTravelTimesGradientVelocity(self):
        """
        Gets the travel times at the receivers in linear gradient velocity model.

        Results
        ttr : matrix
           On successful exit this is the travel times (seconds) from the
           sources to the receivers.  This has dimension [nrec x nsrc]. 
        """ 
        nsrc = self.nsrc
        nrec = self.nrec
        if (nsrc < 1):
            print("No sources")
            return None
        if (nrec < 1): 
            print("No receiver locations set")
            return None
        ttr = ascontiguousarray(zeros(nrec*nsrc), dtype='float64')
        ttrPointer = ttr.ctypes.data_as(POINTER(c_double))
        ierr = c_int(1)
        self.analytic.fteik_analytic_getTravelTimesGradientVel64f(nrec, ttrPointer, ierr)
        if (ierr.value != 0): 
            print("Error getting travel times")
            return None
        if (nsrc > 1):
            ttr = reshape(ttr, [self.nrec, self.nsrc], order='F')
        return ttr

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
        self.analytic.fteik_analytic_setReceivers64f(nrec,
                                                     zrecPointer,
                                                     xrecPointer,
                                                     yrecPointer,
                                                     ierr)
        if (ierr.value != 0): 
            print("Error setting receivers")
            return -1
        self.nrec = nrec
        return 0
