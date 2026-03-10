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

void dif1di(double* f, double kx, double dt, double dx, int nx, int xoff,
	    double alpha)
// Solve 1D diffusion equation using implicit scheme
// f = dependent variable
// kx = diffusion coefficient
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
// alpha: Implicit parameter. 0 for explicit, 0.5 for CN, 1.0 for implicit (default)
{
  int i;
  const double dd=kx*dt/(dx*dx);
  const double add=alpha*dd;
  const double bdd=(1.0-alpha)*dd;
  double a[nx],b[nx],c[nx],d[nx];

  for (i=0;i<nx;i++){
    a[i]=-add;
    b[i]=1.0+2.0*add;
    c[i]=-add;
  }
  for (i=1;i<nx-1;i++){
    d[i]=(1.0-2.0*bdd)*f[i]+bdd*(f[i+1]+f[i-1]);
  }
  trdiag(f,a,b,c,d,nx,xoff);
}

void trdiag(double* s,
	    const double* a, const double* b, const double* c, const double* d,
	    int nx, int xoff)
// Solve a[i]*s[i-1]+b[i]*s[i]+c[i]*s[i+1]=d[i] using Thomas method
{
  int i;
  double bet,gam[nx];

  bet=b[xoff];
  s[xoff]=d[xoff]/bet;
  for (i=xoff+1;i<nx-xoff;i++){
    gam[i]=c[i-1]/bet;
    bet=b[i]-a[i]*gam[i];
    s[i]=(d[i]-a[i]*s[i-1])/bet;
  }
  for (i=nx-xoff-2;i>=xoff;i--) s[i]-=gam[i+1]*s[i+1];
}
