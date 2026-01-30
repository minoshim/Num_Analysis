#include "adv1d.hpp"

void fv3rd(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using 3rd-order finite volume scheme
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  const int sgnv=(v > 0)?1:-1;
  double flux[nx];

  for (i=2;i<nx-1;i++) flux[i]=0.5*(+(1+sgnv)*(-f[i-2]+5*f[i-1]+2*f[i  ])
				    +(1-sgnv)*(-f[i+1]+5*f[i  ]+2*f[i-1]))/6.0; 
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}

void muscl(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using MUSCL scheme
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  const int sgnv=(v > 0)?1:-1;
  double flux[nx];
  double (*msl_func)(double,double,double)=muscl_vl_cal_f;
    
  for (i=2;i<nx-1;i++) flux[i]=0.5*(+(1+sgnv)*msl_func(f[i-2],f[i-1],f[i  ])
				    +(1-sgnv)*msl_func(f[i+1],f[i  ],f[i-1]));
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}
