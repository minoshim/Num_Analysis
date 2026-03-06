import numpy as np
import matplotlib.pyplot as plt

# User-set parameters
## Spatial grid
xoff=2
xmesh=100
nx=xmesh+2*xoff
## Domain size, peak velocity, CFL number
lx=1.0
v=1.0
cfl=0.5
## Simulation time
tmax=0.3
# End of user-set parameters

# Grid width
dx=lx/xmesh
dt=np.abs(cfl*dx/v)
# Variables
x=np.zeros(nx)
f=np.zeros(nx)

# Functions
def init(x,f):                  # Initialize
    nx=len(f)
    i=np.arange(nx)
    x[:]=(i+0.5-xoff)*dx-0.5*lx
    f[:]=np.where((nx//4 < i) & (i < 3*nx//4), v, 0.0)
    # f[:]=v*np.exp(-(x[:]*x[:])/(16*dx*dx))

def bc1d(f,xoff,dnx=0):         # Boundary condition
    nx = len(f)
    if dnx == 0:
        f[:xoff] = f[nx - 2*xoff : nx - xoff]
        f[nx - xoff:] = f[xoff : 2*xoff]
    elif abs(dnx) == 1:
        f[:xoff] = dnx * f[xoff : 2*xoff][::-1]
        f[nx - xoff:] = dnx * f[nx - 2*xoff : nx - xoff][::-1]
        
def bgs_fd(f,dt,dx,xoff=1):
    nx=len(f)
    dtdx=dt/dx
    df=np.zeros_like(f)
    df[1:-1]=np.where(f[1:-1] > 0,(f[1:-1]-f[:-2]),(f[2:]-f[1:-1]))
    f[xoff:nx-xoff]-=dtdx*f[xoff:nx-xoff]*df[xoff:nx-xoff]

def bgs_fv(f,dt,dx,xoff=1):
    nx=len(f)
    dtdx=dt/dx
    flux=np.zeros_like(f)
    flux[1:]=np.where((f[:-1]+f[1:]) > 0,0.5*f[:-1]**2,0.5*f[1:]**2)
    flux[1:]=np.where((f[:-1] < 0) & (f[1:] > 0), 0, flux[1:]) # Correction for expansion wave
    f[xoff:nx-xoff]-=dtdx*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])
    
def main(t,tmax):
    while(t < tmax):
        bc1d(f,xoff,0)
        # bgs_fd(f,dt,dx,xoff)
        bgs_fv(f,dt,dx,xoff)
        t += dt
    return t

if __name__ == "__main__":
    init(x,f)
    t=main(0,tmax)
    print(f"Simulation end at t = {t:.6f}")
