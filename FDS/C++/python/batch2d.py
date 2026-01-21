import numpy as np
import matplotlib.pyplot as plt
from python import plt2d

while True:
    direc=input("Input data directory (Ctrl-D to exit): ")+"/"
    try:
        x=np.loadtxt(direc+"x.dat",dtype=float)
        y=np.loadtxt(direc+"y.dat",dtype=float)
        offs=np.loadtxt(direc+"offs.dat",dtype=int)
        nx=np.size(x)
        ny=np.size(y)
        f=np.loadtxt(direc+"f.dat",dtype=float).reshape((ny,nx))
        rhs=np.loadtxt(direc+"rhs.dat",dtype=float).reshape((ny,nx))
        break
    except:
        print("Error during file load.")

xoff=offs[0]
yoff=offs[1]
dx=x[1]-x[0]
dy=y[1]-y[0]
d2f=np.zeros([ny,nx])
for j in range(1,ny-1):
    for i in range(1,nx-1):
        d2f[j,i]=(+(f[j,i+1]-2*f[j,i]+f[j,i-1])/(dx*dx)
                  +(f[j+1,i]-2*f[j,i]+f[j-1,i])/(dy*dy))

x=x[xoff:nx-xoff]
y=y[yoff:ny-yoff]
f=f[yoff:ny-yoff,xoff:nx-xoff]
rhs=rhs[yoff:ny-yoff,xoff:nx-xoff]
d2f=d2f[yoff:ny-yoff,xoff:nx-xoff]
