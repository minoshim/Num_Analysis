import numpy as np
import matplotlib.pyplot as plt
from python import plt2d

while True:
    direc=input("Input data directory (Ctrl-D to exit): ")+"/"
    try:
        t=np.atleast_1d(np.loadtxt(direc+"t.dat",dtype=float))
        x=np.loadtxt(direc+"x.dat",dtype=float)
        v=np.loadtxt(direc+"v.dat",dtype=float)
        offs=np.loadtxt(direc+"offs.dat",dtype=int)
        nt=np.size(t)
        nx=np.size(x)
        nv=np.size(v)
        f=np.loadtxt(direc+"f.dat",dtype=float).reshape((nt,nv,nx))
        g=np.loadtxt(direc+"g.dat",dtype=float).reshape((nt,nx))
        break
    except:
        print("Error during file load.")

xoff=offs[0]
voff=offs[1]
dx=x[1]-x[0]
dv=v[1]-v[0]

x=x[xoff:nx-xoff]
v=v[voff:nv-voff]
f=f[:,voff:nv-voff,xoff:nx-xoff]
g=g[:,xoff:nx-xoff]

# Potential
p=np.zeros((nt,nx-2*xoff))
for i in range(nx-2*xoff-1):
    p[:,i+1]=p[:,i]-g[:,i+1]*dx
p=p-p.mean(axis=1,keepdims=True)

# Total number, energy, and entropy
num=f.sum(axis=(1,2))*dx*dv
ene=0.5*((f*v[None,:,None]*v[None,:,None]).sum(axis=(1,2))+(g*g).sum(axis=1))*dx*dv
ent=(f*(1-f)).sum(axis=(1,2))*dx*dv
