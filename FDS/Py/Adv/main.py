import numpy as np
import matplotlib.pyplot as plt

xoff=2
xmesh=100
nx=xmesh+2*xoff
lx=1.0
dx=lx/xmesh
v=1.0
cfl=0.5
dt=np.abs(cfl*dx/v)

if __name__ == "__main__":
    print(0)
    
