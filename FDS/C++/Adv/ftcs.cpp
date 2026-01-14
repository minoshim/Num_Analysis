#include "ftcs.hpp"

void ftcs(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using FTCS scheme
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  double df[nx];

  for (i=1;i<nx-1;i++) df[i]=0.5*(f[i+1]-f[i-1]);
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*df[i];
}
