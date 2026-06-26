import numpy as np
import matplotlib.pyplot as plt

# User-set parameters
## Spatial grid
xoff=2
xmesh=100
nx=xmesh+2*xoff
## Domain size, advection velocity, CFL number
lx=1.0
v=1.0
cfl=0.5
## Simulation time
tmax=0.15
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
    f[:]=np.where((nx//4 < i) & (i < 3*nx//4), 1.0, 0.0)
    # f[:]=np.exp(-(x[:]*x[:])/(16*dx*dx))

def bc1d(f,xoff,dnx=0):         # Boundary condition
    nx = len(f)
    if dnx == 0:
        f[:xoff] = f[nx - 2*xoff : nx - xoff]
        f[nx - xoff:] = f[xoff : 2*xoff]
    elif abs(dnx) == 1:
        f[:xoff] = dnx * f[xoff : 2*xoff][::-1]
        f[nx - xoff:] = dnx * f[nx - 2*xoff : nx - xoff][::-1]

def minmod(a,b):
    return np.where(a*b < 0, 0, np.where(np.abs(a) < np.abs(b), a, b))

def median3(a,b,c):
    return np.maximum(np.minimum(b,c),np.minimum(a,np.maximum(b,c)))

def musclmm(f,v,dt,dx,xoff=2):
    nx=len(f)
    nu=v*dt/dx
    sgnv=np.sign(v)
    flux=np.zeros_like(f)
    slopel=minmod(f[1:-2]-f[0:-3],f[2:-1]-f[1:-2])
    sloper=minmod(f[2:-1]-f[1:-2],f[3:  ]-f[2:-1])
    flux[2:-1]=0.5*(+(1+sgnv)*(f[1:-2]+0.5*slopel)+(1-sgnv)*(f[2:-1]-0.5*sloper))
    f[xoff:nx-xoff]-=nu*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])
    
def fv3rd(f,v,dt,dx,xoff=2):       # 3rd-order finite-volume scheme
    nx=len(f)
    nu=v*dt/dx
    sgnv=np.sign(v)
    flux=np.zeros_like(f)
    flux[2:-1]=0.5*(+(1+sgnv)*(-f[0:-3]+5*f[1:-2]+2*f[2:-1])
                    +(1-sgnv)*(-f[3:  ]+5*f[2:-1]+2*f[1:-2]))/6.0
    f[xoff:nx-xoff]-=nu*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])    

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

def pfc_flux(f0,a1,a2,xi,sgnv):
    c0=f0-a2/12.0
    if sgnv > 0:
        xm=0.5*(1.0-xi)
        xl=0.5-xi
        return (c0+0.5*a1+0.25*a2
                +4.0*(c0+a1*xm+a2*xm*xm)
                +c0+a1*xl+a2*xl*xl)/6.0
    xm=-0.5*(1.0-xi)
    xr=-0.5+xi
    return (c0-0.5*a1+0.25*a2
            +4.0*(c0+a1*xm+a2*xm*xm)
            +c0+a1*xr+a2*xr*xr)/6.0

def pfc_coef(fm,f0,fp):
    fmax=np.maximum(np.maximum(fm,fp),np.minimum(2.0*f0-fm,2.0*f0-fp))
    fmin=np.maximum(0.0,np.minimum(np.minimum(fm,fp),np.maximum(2.0*f0-fm,2.0*f0-fp)))
    alpha=1.0/3.0
    ap=fp-f0
    am=f0-fm
    fpmin=3.0*np.maximum(2.0*(f0-fmax),fmin-f0)
    fpmax=3.0*np.minimum(2.0*(f0-fmin),fmax-f0)
    fmmin=-fpmax
    fmmax=-fpmin
    ap=median3(ap,alpha*fpmin,alpha*fpmax)
    am=median3(am,alpha*fmmin,alpha*fmmax)
    a1=0.5*(ap+am)
    a2=0.5*(ap-am)
    return a1,a2

def pfc(f,v,dt,dx,xoff=2):       # positive and flux conservative scheme
    nx=len(f)
    nu=v*dt/dx
    anu=np.abs(nu)
    sgnv=np.sign(v)
    flux=np.zeros_like(f)
    if anu == 0.0:
        return
    if sgnv > 0:
        fm=f[0:-3]
        f0=f[1:-2]
        fp=f[2:-1]
    else:
        fm=f[1:-2]
        f0=f[2:-1]
        fp=f[3:  ]
    a1,a2=pfc_coef(fm,f0,fp)
    flux[2:-1]=pfc_flux(f0,a1,a2,anu,sgnv)
    f[xoff:nx-xoff]-=nu*(flux[xoff+1:nx-xoff+1]-flux[xoff:nx-xoff])
    
def main(t,tmax):
    fcpy=np.zeros_like(f)
    while(t < tmax):
        bc1d(f,xoff,0)
        if (1):                 # 1 for SL scheme, 0 for FV-RK2 scheme
            csl3rd(f,v,dt,dx,xoff)
            # pfc(f,v,dt,dx,xoff)
        else :
            fcpy=f.copy()
            musclmm(f,v,dt,dx,xoff)
            bc1d(f,xoff,0)
            musclmm(f,v,dt,dx,xoff)
            f[xoff:nx-xoff]=0.5*(f[xoff:nx-xoff]+fcpy[xoff:nx-xoff])
        t += dt
    return t
    
if __name__ == "__main__":
    init(x,f)
    t=main(0,tmax)
    print(f"Simulation end at t = {t:.6f}")
