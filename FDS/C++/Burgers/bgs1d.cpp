#include "bgs1d.hpp"

void bgs_fd(double* f, double dt, double dx, int nx, int xoff)
// Solve 1D Burgers equation using non-conservative finite-difference scheme
// f = dependent variable
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double dtdx=dt/dx;
  int sgnv;
  double df[nx];

  for (i=1;i<nx-1;i++){
    sgnv=(f[i] > 0)?1:-1;
    df[i]=0.5*((1+sgnv)*(f[i]-f[i-1])+(1-sgnv)*(f[i+1]-f[i]));    
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=dtdx*f[i]*df[i];
}

void bgs_fv(double* f, double dt, double dx, int nx, int xoff)
// Solve 1D Burgers equation using conservative finite-volume scheme
// f = dependent variable
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double dtdx=dt/dx;
  int sgnv;
  double flux[nx];

  for (i=1;i<nx-1;i++){
    sgnv=((f[i-1]+f[i]) > 0)?1:-1;
    flux[i]=0.25*((1+sgnv)*(f[i-1]*f[i-1])+(1-sgnv)*(f[i]*f[i]));
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=dtdx*(flux[i+1]-flux[i]);
}
