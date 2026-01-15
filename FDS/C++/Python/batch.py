import numpy as np
import matplotlib.pyplot as plt

while True:
    direc=input("Input data directory (Ctrl-D to exit): ")+"/"
    try:
        t=np.atleast_1d(np.loadtxt(direc+"t.dat",dtype=float))
        x=np.loadtxt(direc+"x.dat",dtype=float)
        xoff=int(np.loadtxt(direc+"xoff.dat",dtype=float))
        nt=np.size(t)
        nx=np.size(x)
        f=np.loadtxt(direc+"f.dat",dtype=float).reshape((nt,nx))
        break
    except:
        print("Error during file load.")

x=x[xoff:nx-xoff]
f=f[:,xoff:nx-xoff]
