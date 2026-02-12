#include "adv1d.hpp"

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

void upwd(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using upwind scheme
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  const int sgnv=(v > 0)?1:-1;
  double df[nx];

  for (i=1;i<nx-1;i++) df[i]=0.5*((1+sgnv)*(f[i]-f[i-1])+(1-sgnv)*(f[i+1]-f[i]));
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*df[i];
}

void lawe(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using Lax-Wendroff scheme
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  double df[nx];
  
  for (i=1;i<nx-1;i++) df[i]=0.5*((f[i+1]-f[i-1])-nu*(f[i+1]-2*f[i]+f[i-1]));
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*df[i];
}

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

void sl3rd(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using 3rd-order semi-Lagrange scheme (CFL<1 only)
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  const int sgnv=(v > 0)?1:-1;
  const double anu=sgnv*nu;
  double fcpy[nx];

  for (i=0;i<nx;i++) fcpy[i]=f[i];
  for (i=2;i<nx-2;i++){
    double c0=fcpy[i];
    double c1=(-fcpy[i-2*sgnv]+6*fcpy[i-sgnv]-3*fcpy[i]-2*fcpy[i+sgnv])/6;
    double c2=(               +  fcpy[i-sgnv]-2*fcpy[i]+  fcpy[i+sgnv])*0.5;
    double c3=(+fcpy[i-2*sgnv]-3*fcpy[i-sgnv]+3*fcpy[i]-  fcpy[i+sgnv])/6;
    f[i]=c0+anu*(c1+anu*(c2+anu*c3));
  }
}
