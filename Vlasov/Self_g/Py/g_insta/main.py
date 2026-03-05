import numpy as np
import matplotlib.pyplot as plt

# User-set parameters
## Spatial grid
xoff=2
voff=2
xmesh=64
vmesh=64
nx=xmesh+2*xoff
nv=vmesh+2*voff
## Domain size and CFL number
lx=1.0                          # X domain, [-lx/2,lx/2]
lv=1.0                          # V domain, [-lv,lv]
cfl=0.2
## Simulation time
tmax=3.0
## Perturbation
amp=0.01                        # Perturbation amplitude
k0=2.0*np.pi/lx                 # Fundamental wavenumber
kk=k0*2.0                       # Perturbation wavenumber
kj=kk/0.5                       # Jeans wavenumber
# End of user-set parameters

# Grid width and thermal velocity
dx=lx/xmesh
dv=2*lv/vmesh
dt=np.abs(cfl*dx/lv)
vs=np.sqrt(4*np.pi)/kj
# Variables
x=np.zeros(nx)
v=np.zeros(nv)
f=np.zeros((nv,nx))
g=np.zeros(nx)
drho=np.zeros(nx)

# Functions
def gaussian(x,a,b):
    return np.exp(-0.5*((x-a)/b)**2) / (b*np.sqrt(2*np.pi))

def init(x,v,f):                  # Initialize
    nv,nx=f.shape
    i=np.arange(nx)
    j=np.arange(nv)
    x[:]=(i+0.5-xoff)*dx-0.5*lx
    v[:]=(j+0.5-voff)*dv-lv
    fv=gaussian(v,0,vs)
    fx=1.0+amp*np.cos(kk*x)
    f[:,:]=fv[:,None]*fx[None,:]

def bc2d(f,xoff,yoff,dnx=0,dny=0): # Boundary condition
    # f : shape (ny,nx)
    # dnx, dny: 0 (periodic), +-1 (free of fix)
    ny, nx = f.shape
    # ---- x direction ----
    if dnx == 0:
        f[:, :xoff]      = f[:, nx-2*xoff:nx-xoff]
        f[:, nx-xoff:]   = f[:, xoff:2*xoff]
    elif abs(dnx) == 1:
        f[:, :xoff] = dnx * f[:, xoff:2*xoff][:, ::-1] # Left
        f[:, nx-xoff:] = dnx * f[:, nx-2*xoff:nx-xoff][:, ::-1] # Right
    # ---- y direction ----
    if dny == 0:
        f[:yoff, :]      = f[ny-2*yoff:ny-yoff, :]
        f[ny-yoff:, :]   = f[yoff:2*yoff, :]
    elif abs(dny) == 1:
        f[:yoff, :] = dny * f[yoff:2*yoff, :][::-1, :] # Bottom
        f[ny-yoff:, :] = dny * f[ny-2*yoff:ny-yoff, :][::-1, :] # Top

def csl3rd(f,v,dt,dx,xoff=2):       # 3rd-order conservative semi-Lagrangian scheme
    nx=len(f)
    nu=v*dt/dx
    sgnv=np.sign(v)
    flux=np.zeros_like(f)
    c0l=(-f[0:-3]+5*f[1:-2]+2*f[2:-1])/6.0
    c0r=(-f[3:  ]+5*f[2:-1]+2*f[1:-2])/6.0
    c1 =f[2:-1]-f[1:-2]
    c2l=(f[0:-3]-2*f[1:-2]+f[2:-1])*0.5
    c2r=(f[3:  ]-2*f[2:-1]+f[1:-2])*0.5
    ftl=c0l+nu*(-c1*0.5+c2l*nu/3.0)
    ftr=c0r+nu*(-c1*0.5+c2r*nu/3.0)
    flux[2:-1]=0.5*((1+sgnv)*ftl+(1-sgnv)*ftr)
    f[xoff:nx-xoff]-=nu*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])    

def pushx(f,v,dt,dx,xoff=2,voff=2):
    bc2d(f,xoff,voff,0,999) # Periodic in X, nothing to do in V
    for fi,vi in zip(f[voff:-voff,:],v[voff:-voff]):
        csl3rd(fi,vi,dt,dx,xoff)

def pushv(f,g,dt,dv,xoff=2,voff=2):
    nv,nx=f.shape
    for i in range(xoff,nx-xoff):
        ftmp = f[:, i].copy()
        csl3rd(ftmp, 0.5*(g[i]+g[i+1]), dt, dv, voff)
        f[voff:nv-voff, i] = ftmp[voff:nv-voff]
    
def main(t,tmax):
    while(t < tmax):
        # Advection in X (half)
        pushx(f,v,0.5*dt,dx,xoff,voff)

        # Calculate density fluctuation
        drho[xoff:-xoff]=np.sum(f[voff:-voff,xoff:-xoff],axis=0)*dv
        dmean=np.sum(drho[xoff:-xoff])*dx/lx
        drho[xoff:-xoff]-=dmean

        # Calculate gravity field
        g[xoff]=0.0
        for i in range(xoff+1,nx-xoff+1):
            g[i]=g[i-1]-4.0*np.pi*drho[i-1]*dx
        gmean=np.sum(0.5*(g[xoff:nx-xoff]+g[xoff+1:nx-xoff+1]))*dx/lx
        g[xoff:nx-xoff+1]-=gmean
        
        # Advection in V (full)
        pushv(f,g,    dt,dv,xoff,voff)
        
        # Advection in X (half)
        pushx(f,v,0.5*dt,dx,xoff,voff)
            
        t+=dt
    return t
    
if __name__ == "__main__":
    init(x,v,f)
    t=main(0,tmax)
    print(f"Simulation end at t = {t:.6f}")
