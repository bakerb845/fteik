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
from fteikAnalytic import fteikAnalytic
from fteik2d import fteik2d
from fteik3d import fteik3d
import os


class fteik:
    def __init__(self,
                 fteik_library='/usr/local/lib/libfteik_shared.so'):
        self.fteik2d = fteik2d(fteik_library)
        self.fteik3d = fteik3d(fteik_library)

    def __enter__(self):
        return self

    def __exit__(self):
        self.fteik2d.free()
        self.fteik3d.free()
        return

    def free(self):
        """
        Frees the internal Fortran module.
        """
        self.fteik2d.free()
        self.fteik3d.free()
        return

if __name__ == "__main__":
    from matplotlib import pyplot as plt
    nx = 101
    nz = 101
    dx = 10.0
    dz = 10.0
    nrec = 11
    fteikLibrary = '/home/bakerb25/C/fteik/lib/libfteik_shared.so'
    x = linspace(0.0, (nx-1)*dx, nx)
    z = linspace(0.0, (nz-1)*dz, nz)
    xsrc = (nx - 1)*dx/4.0
    zsrc = (nz - 1)*dz/4.0
    xrec = linspace((nx-1)*dx*0.2, (nx-1)*dx*0.8, nrec)
    zrec = zeros(nrec)
    vel = zeros([nz,nx]) + 5.e3 #zeros([nz-1,nx-1]) + 5.e3
    fteik2d = fteik2d(fteikLibrary)
    fteik2d.initialize(nx, nz, dx, dz, verbose=3)
    fteik2d.setVelocityModel(vel) 
    fteik2d.setSources(xsrc, zsrc)
    fteik2d.setReceivers(xrec, zrec)
    fteik2d.solveLSM() 
    ttimes = fteik2d.getTravelTimeField()
    ttr = fteik2d.getTravelTimes()
    print(ttr)
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
    fteik3d = fteik3d(fteikLibrary)
    fteik3d.initialize(nx, ny, nz, dx, dy, dz, verbose=3) 
    fteik3d.setVelocityModel(vel)
    fteik3d.setSources(xsrc, ysrc, zsrc)
    fteik3d.solveLSM()
    ttimes = fteik3d.getTravelTimeField()
    fteik3d.free()

    analytic = fteikAnalytic(fteikLibrary)
    analytic.initialize(nx, ny, nz, dx, dy, dz, verbose=3)
    #analytic.setVelocityModel(vel
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
