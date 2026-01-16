#include "dif1d.hpp"

void dif1de(double* f, double kx, double dt, double dx, int nx, int xoff)
// Solve 1D diffusion equation using explicit scheme
// f = dependent variable
// kx = diffusion coefficient
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double dd=kx*dt/(dx*dx);
  double df[nx];

  for (i=1;i<nx;i++) df[i]=f[i]-f[i-1];
  for (i=xoff;i<nx-xoff;i++) f[i]+=dd*(df[i+1]-df[i]);
}
