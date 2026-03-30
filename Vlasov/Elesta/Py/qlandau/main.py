import numpy as np
import matplotlib.pyplot as plt

# Quantum Vlasov simulation

# User-set parameters
## Spatial grid
xoff=2
voff=2
xmesh=128
vmesh=128
nx=xmesh+2*xoff
nv=vmesh+2*voff
## Domain size and CFL number
lx=4.0*np.pi                          # X domain, [-lx/2,lx/2]
lv=5.0                          # V domain, [-lv,lv]
cfl=0.2
## Simulation time
tmax=16.0
# Normalized charge
q=-1.0                          # -1 for electrons
## Perturbation
amp=0.05                        # Perturbation amplitude
kk=0.5                       # Perturbation wavenumber
vd=1.5                       # Drift velocity for two-stream instability
vs=0.2                       # Thermal velocity for two-stream instability
## Quantum parameters
hbar=1e0
hbar2=hbar*hbar/24.0
# End of user-set parameters

# Grid width
dx=lx/xmesh
dv=2*lv/vmesh
dt=np.abs(cfl*dx/lv)
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
    # fv=gaussian(v,0,1.0)
    fv=0.5*(gaussian(v,+vd,vs)+gaussian(v,-vd,vs)) # Two-stream instability
    fx=1.0+amp*np.cos(kk*x)
    f[:,:]=fv[:,None]*fx[None,:]

def bc1d(f,xoff,dnx=0):         # Boundary condition
    nx = len(f)
    if dnx == 0:
        f[:xoff] = f[nx - 2*xoff : nx - xoff]
        f[nx - xoff:] = f[xoff : 2*xoff]
    elif abs(dnx) == 1:
        f[:xoff] = dnx * f[xoff : 2*xoff][::-1]
        f[nx - xoff:] = dnx * f[nx - 2*xoff : nx - xoff][::-1]

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

def minmod(a,b):
    return np.where(a*b < 0, 0, np.where(np.abs(a) < np.abs(b), a, b))
def mm3(a,b,c):
    return minmod(minmod(a,b),c)
def cslmsl(f,v,dt,dx,xoff=2):       # Conservative SL-MUSCL scheme
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
    slopel=mm3(f[1:-2]-f[0:-3],f[2:-1]-f[1:-2],ftl-f[1:-2])
    sloper=mm3(f[2:-1]-f[1:-2],f[3:  ]-f[2:-1],f[2:-1]-ftr)
    flux[2:-1]=0.5*((1+sgnv)*(f[1:-2]+slopel)+(1-sgnv)*(f[2:-1]-sloper))
    rhs=-nu*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])
    f[xoff:nx-xoff]+=rhs
    return rhs
        
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
    rhs=-nu*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])    
    f[xoff:nx-xoff]+=rhs
    return rhs

def dadv1d(f,v,w,dt,dx,xoff=2,alpha=1.0):
    # Solve df/dt+v*df/dx+w*d3f/dx3=0
    # Fully implicit scheme
    from scipy.linalg import solve_banded
    nx=len(f)
    vv =v*dt/dx                  # Advection
    avv=alpha*vv
    dd =0.5*np.abs(vv)           # Diffusion (numerical)
    add=alpha*dd
    ww =w*dt/(dx*dx*dx)          # Dispersion
    aww=alpha*ww
    ss =0.5*np.abs(ww)          # Hyper diffusion (numerical)
    ass=alpha*ss
    c1=np.array([-1./12.,+2./3.,+0.0,  -2./3.,+1./12.])
    c2=np.array([-1./6. ,+2./3.,-1.0,  +2./3., -1./6.]) # These coefficients give 3rd-order upwind scheme for advection
    c3=np.array([+0.5   ,-1.0  ,+0.0,    +1.0,   -0.5])
    c4=np.array([+1.0   ,-4.0  ,+6.0,    -4.0,   +1.0])
    ab=np.zeros((5,nx)) #banded matrix
    ab[0, 2:]=c1[0]*avv-c2[0]*add+c3[0]*aww+c4[0]*ass #i+2
    ab[1, 1:]=c1[1]*avv-c2[1]*add+c3[1]*aww+c4[1]*ass #i+1
    ab[2,  :]=1.0      -c2[2]*add          +c4[2]*ass #i
    ab[3,:-1]=c1[3]*avv-c2[3]*add+c3[3]*aww+c4[3]*ass #i-1
    ab[4,:-2]=c1[4]*avv-c2[4]*add+c3[4]*aww+c4[4]*ass #i-2
    rr=np.zeros(nx)     #RHS
    rr[2:-2]=(-vv*(c1[0]*f[4:]+c1[1]*f[3:-1]+c1[2]*f[2:-2]+c1[3]*f[1:-3]+c1[4]*f[0:-4])
              +dd*(c2[0]*f[4:]+c2[1]*f[3:-1]+c2[2]*f[2:-2]+c2[3]*f[1:-3]+c2[4]*f[0:-4])
              -ww*(c3[0]*f[4:]+c3[1]*f[3:-1]+c3[2]*f[2:-2]+c3[3]*f[1:-3]+c3[4]*f[0:-4])
              -ss*(c4[0]*f[4:]+c4[1]*f[3:-1]+c4[2]*f[2:-2]+c4[3]*f[1:-3]+c4[4]*f[0:-4]))
    f[xoff:-xoff]+=solve_banded((2,2),ab[:,xoff:-xoff],rr[xoff:-xoff])

def dcsl1d(f,v,w,dt,dx,xoff=2,alpha=1.0):
    # Solve df/dt+v*df/dx+w*d3f/dx3=0
    # CSL scheme for advection, implicit scheme for dispersion
    from scipy.linalg import solve_banded
    nx=len(f)
    ww =w*dt/(dx*dx*dx)          # Dispersion
    aww=alpha*ww
    ss =0.5*np.abs(ww)          # Hyper diffusion (numerical)
    ass=alpha*ss
    c3=np.array([+0.5   ,-1.0  ,+0.0,    +1.0,   -0.5])
    c4=np.array([+1.0   ,-4.0  ,+6.0,    -4.0,   +1.0])
    ab=np.zeros((5,nx)) #banded matrix
    ab[0, 2:]=c3[0]*aww+c4[0]*ass #i+2
    ab[1, 1:]=c3[1]*aww+c4[1]*ass #i+1
    ab[2,  :]=1.0      +c4[2]*ass #i
    ab[3,:-1]=c3[3]*aww+c4[3]*ass #i-1
    ab[4,:-2]=c3[4]*aww+c4[4]*ass #i-2
    rr=np.zeros(nx)     #RHS
    fcpy=f.copy()       # Dummy
    rr[2:-2]=(cslmsl(fcpy,v,dt,dx,2)
              -ww*(c3[0]*f[4:]+c3[1]*f[3:-1]+c3[2]*f[2:-2]+c3[3]*f[1:-3]+c3[4]*f[0:-4])
              -ss*(c4[0]*f[4:]+c4[1]*f[3:-1]+c4[2]*f[2:-2]+c4[3]*f[1:-3]+c4[4]*f[0:-4]))
    f[xoff:-xoff]+=solve_banded((2,2),ab[:,xoff:-xoff],rr[xoff:-xoff])
    
def pushx(f,v,dt,dx,xoff=2,voff=2):
    bc2d(f,xoff,voff,0,999) # Periodic in X, nothing to do in V
    for fi,vi in zip(f[voff:-voff,:],v[voff:-voff]):
        cslmsl(fi,vi,dt,dx,xoff)

def pushv(f,g,dt,dv,dx,xoff=2,voff=2):
    bc1d(g,xoff,0)
    nv,nx=f.shape
    for i in range(xoff,nx-xoff):
        ftmp = f[:, i].copy()
        d2gdx=0.5*(g[i+2]-g[i+1]-g[i]+g[i-1])/(dx*dx)
        dcsl1d(ftmp, 0.5*(g[i]+g[i+1]), hbar2*d2gdx, dt, dv, voff, 0.505)
        # dadv1d(ftmp, 0.5*(g[i]+g[i+1]), hbar2*d2gdx, dt, dv, voff, 0.505)
        # cslmsl(ftmp, 0.5*(g[i]+g[i+1]), dt, dv, voff)
        f[voff:nv-voff, i] = ftmp[voff:nv-voff]
    
def main(t,tmax):
    while(t < tmax):
        # Advection in X (half)
        pushx(f,v,0.5*dt,dx,xoff,voff)

        # Calculate density fluctuation
        drho[xoff:-xoff]=np.sum(f[voff:-voff,xoff:-xoff],axis=0)*dv
        dmean=np.sum(drho[xoff:-xoff])*dx/lx
        drho[xoff:-xoff]-=dmean

        # Calculate electric field
        g[xoff]=0.0
        for i in range(xoff+1,nx-xoff+1):
            g[i]=g[i-1]+q*drho[i-1]*dx
        gmean=np.sum(0.5*(g[xoff:nx-xoff]+g[xoff+1:nx-xoff+1]))*dx/lx
        g[xoff:nx-xoff+1]-=gmean

        # Advection in V (full)
        pushv(f,q*g,  dt,dv,dx,xoff,voff)
        
        # Advection in X (half)
        pushx(f,v,0.5*dt,dx,xoff,voff)
            
        t+=dt
    return t
    
if __name__ == "__main__":
    init(x,v,f)
    t=main(0,tmax)
    print(f"Simulation end at t = {t:.6f}")
    plt.imshow(f[voff:-voff,xoff:-xoff],cmap='jet',origin='lower',aspect='auto',extent=[-lx/2,lx/2,-lv,lv])
    plt.show(block=False)
