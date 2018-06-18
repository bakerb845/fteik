#!/usr/bin/env python3
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
from numpy import sqrt
from numpy import ascontiguousarray 
from numpy import linspace
from fteikAnalytic import fteikAnalytic
from fteik2d import fteik2d
from fteik3d import fteik3d
import os


class fteik:
    def __init__(self,
                 fteik_path=os.environ['LD_LIBRARY_PATH'].split(os.pathsep),
                 fteik_library='libfteik_shared.so'):
        #         fteik_library='/usr/local/lib/libfteik_shared.so'):
        self.fteik2d = fteik2d(fteik_path, fteik_library)
        self.fteik3d = fteik3d(fteik_path, fteik_library)

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
    ny = 1
    nz = 101
    dx = 10.0
    dy = 1.0
    dz = 10.0
    nrec = 11
    vconst = 5.e3
    fteikLibrary = '/home/bakerb25/C/fteik/lib/libfteik_shared.so'
    x = linspace(0.0, (nx-1)*dx, nx)
    z = linspace(0.0, (nz-1)*dz, nz)
    xsrc = (nx - 1)*dx/4.0
    ysrc = 0.0
    zsrc = (nz - 1)*dz/4.0
    xrec = linspace((nx-1)*dx*0.2, (nx-1)*dx*0.8, nrec)
    yrec = zeros(nrec)
    zrec = zeros(nrec)
    vmin = vconst
    vmax = vconst + 1000.0
    print(xsrc, zsrc)
    vel = zeros([nz,nx]) + vconst
    fteik2d = fteik2d()#fteikLibrary)
    fteik2d.initialize(nx, nz, dx, dz, verbose=3)
    fteik2d.setVelocityModel(vel) 
    fteik2d.setSources(xsrc, zsrc)
    fteik2d.setReceivers(xrec, zrec)
    fteik2d.solveLSM() 
    ttimesConst = fteik2d.getTravelTimeField()
    ttrConst = fteik2d.getTravelTimes()
    #print(ttrConst)
    fteik2d.free()
    # Get the analytic times 
    analytic = fteikAnalytic()#fteikLibrary)
    analytic.initialize(nx, ny, nz, dx, dy, dz, verbose=3)
    analytic.setSources(xsrc, ysrc, zsrc)
    analytic.setReceivers(xrec, yrec, zrec)

    analytic.setConstantVelocity(vconst)
    analytic.solveConstantVelocity()
    ttConstantAnalytic = analytic.getTravelTimeField() 
    ttrConstantAnalytic = analytic.getTravelTimesConstantVelocity()

    #analytic.setLinearVelocityGradient(vmin, vmax)
    #analytic.solveLinearVelocityGradient()
    #ttGradAnalytic = analytic.getTravelTimeField()
    #print(ttGradAnalytic.T)
    #ttrGradAnalytic = analytic.getTravelTimesGradientVelocity()
    analytic.free()
    tres = ttConstantAnalytic - ttimesConst
    tresL1 = abs(tres).sum()/tres.size
    tresL8 = abs(tres).max()
    tresL2 = sqrt(pow(tres, 2).sum())/tres.size

    print("Gradient L1/L2/Linf residual:", tresL1, tresL2, tresL8)
    print("Constant travel times\n", ttrConst)
    print("Reference constant travel times\n", ttrConstantAnalytic)

    #print("Gradient travel times\n", ttrGrad)
    #print("Reference gradient travel times\n", ttrGradAnalytic)
    #stop

    # plot it
    """
    ttimes = ttimesConst
    print(ttimes.shape)
    plt.contourf(x, z, ttimes, 25, cmap=plt.cm.viridis)
    plt.ylim((nz - 1)*dz, 0)
    plt.title('Travel Time Field (s)')
    plt.xlabel('Offset (m)')
    plt.ylabel('Depth (m)')
    plt.colorbar()
    plt.show()
    stop
    """

    nx = 90
    ny = 90
    nz = 90
    dx = 15.0
    dy = 15.0
    dz = 15.0
    xsrc = (nx - 1)*dx/4.0
    ysrc = (ny - 1)*dy/2.0
    zsrc = (nz - 1)*dz/4.0
    xrec = zeros(nrec) + (nx - 1)*dx/2.0
    yrec = zeros(nrec) + (ny - 1)*dy/2.0
    zrec = linspace((nz-1)*dz*0.2, (nz-1)*dz*0.8, nrec)
    vmin = 4000.0
    vmax = 6000.0 
    vel = zeros([nz,ny,nx]) + vconst
    velGrad = zeros([nz,ny,nx]) + vconst
    # linear velocity gradient model 
    for iz in range(nz):
        velGrad[iz,:,:] = vmin + (vmax - vmin)/((nz-1)*dz)*(iz*dz)
        #print(vmin + (vmax - vmin)/((nz-1)*dz)*(iz*dz))
    fteik3d = fteik3d()#fteikLibrary)
    fteik3d.initialize(nx, ny, nz, dx, dy, dz, verbose=3) 
    fteik3d.setSources(xsrc, ysrc, zsrc)
    fteik3d.setReceivers(xrec, yrec, zrec)
    # Solve the constant velocity problem
    fteik3d.setVelocityModel(vel)
    fteik3d.solveLSM()
    ttimesConst = fteik3d.getTravelTimeField() 
    ttrConst = fteik3d.getTravelTimes()
    # Solve the linear gradient problem
    fteik3d.setVelocityModel(velGrad)
    fteik3d.solveLSM()
    ttimesGrad = fteik3d.getTravelTimeField()
    ttrGrad = fteik3d.getTravelTimes()
    fteik3d.free()

    analytic = fteikAnalytic()#fteikLibrary)
    analytic.initialize(nx, ny, nz, dx, dy, dz, verbose=3)
    analytic.setSources(xsrc, ysrc, zsrc)
    analytic.setReceivers(xrec, yrec, zrec)
    analytic.setConstantVelocity(vconst)
    analytic.setLinearVelocityGradient(vmin, vmax)
    analytic.solveConstantVelocity()
    ttConstantAnalytic = analytic.getTravelTimeField() 
    ttrConstantAnalytic = analytic.getTravelTimesConstantVelocity()
    analytic.solveLinearVelocityGradient()
    ttGradAnalytic = analytic.getTravelTimeField()
    ttrGradAnalytic = analytic.getTravelTimesGradientVelocity()
    analytic.free()
    tres = ttConstantAnalytic - ttimesConst
    tresL1 = abs(tres).sum()/tres.size
    tresL8 = abs(tres).max()        
    tresL2 = sqrt(pow(tres, 2).sum())/tres.size 
    print("Constant L1/L2/Linf residual:", tresL1, tresL2, tresL8)
    tres = ttGradAnalytic - ttimesGrad
    tresL1 = abs(tres).sum()/tres.size
    tresL8 = abs(tres).max()
    tresL2 = sqrt(pow(tres, 2).sum())/tres.size

    print("Gradient L1/L2/Linf residual:", tresL1, tresL2, tresL8)
    print("Constant travel times\n", ttrConst)
    print("Reference constant travel times\n", ttrConstantAnalytic)
    print("Gradient travel times\n", ttrGrad)
    print("Reference gradient travel times\n", ttrGradAnalytic)

    #, tres.max(), ttGradAnalytic.max(), ttimesGrad.max())


    #"""
    ttimes = ttimesGrad#ConstantAnalytic
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
    #"""
