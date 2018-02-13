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
        lib.fteik_analytic_setVelocityGradient64f.argtypes = (c_double,      #v0
                                                              c_double,      #v1
                                                              POINTER(c_int) #ierr
                                                             )
        lib.fteik_analytic_setVelocityConstant64f.argtypes = (c_double,      #v0 
                                                              POINTER(c_int) #ierr
                                                             )
        lib.fteik_analytic_getTravelTimeField64f.argtypes = (c_int,             #ngrd
                                                             c_int,             #order
                                                             POINTER(c_double), #ttimes
                                                             POINTER(c_int)     #ierr
                                                            )
        lib.fteik_analytic_free.argtypes = None

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


