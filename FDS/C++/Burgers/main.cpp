#include "cnst.hpp"
#include "routine.hpp"
#include "bgs1d.hpp"

int main(void)
{
  int i,n;
  const int xoff=XOFF;
  const int nx=XMESH+2*xoff;
  const int nrec=NREC;
  const int nmax=NMAX;
  const double lx=1.0;		// Domain size
  const double dx=lx/XMESH;	// Grid spacing
  const double v=1.0;		// Peak Velcity
  const double dt=fabs(CFL*dx/v); // Time step
  double t=0.0;
  double x[nx],f[nx];

  // Initialize
  for (i=0;i<nx;i++){
    x[i]=(i+0.5-xoff)*dx-0.5*lx;
    f[i]=(i > nx/4 && i < 3*nx/4)?v:0; // Square wave
    // f[i]=v*exp(-(x[i]*x[i])/(64*dx*dx));	   // Gaussian
  }

  // Output
  fout(&xoff,"xoff.dat",1,0);
  fout( x,"x.dat",nx,0);
  fout( f,"f.dat",nx,0);
  fout(&t,"t.dat",1, 0);
  
  // Integration
  while(n++ < nmax){
    bc1d(f,nx,xoff,0);
    // bgs_fd(f,dt,dx,nx,xoff);
    bgs_fv(f,dt,dx,nx,xoff);
    t+=dt;

    // Output
    if (n % nrec == 0){
      fout( f,"f.dat",nx,1);
      fout(&t,"t.dat",1, 1);
    }
  }
  
  return 0;
}
