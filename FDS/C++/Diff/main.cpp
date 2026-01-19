#include "cnst.hpp"
#include "routine.hpp"

int main(void)
{
  int i,n;
  const int xoff=XOFF;
  const int nx=XMESH+2*xoff;
  const int nrec=NREC;
  const int nmax=NMAX;
  const double lx=1.0;		// Domain size
  const double dx=lx/XMESH;	// Grid spacing
  const double kx=1.0;		// Diffusion coefficient
  const double dt=fabs(CFL*dx*dx/kx); // Time step
  double t=0.0;
  double *x,*f;
  x=new double[nx]();
  f=new double[nx]();

  // Initialize
  for (i=0;i<nx;i++){
    x[i]=(i+0.5-xoff)*dx-0.5*lx;
    f[i]=exp(-(x[i]*x[i])/(16*dx*dx));	   // Gaussian
  }

  // Output
  fout(&xoff,"xoff.dat",1,0);
  fout( x,"x.dat",nx,0);
  fout( f,"f.dat",nx,0);
  fout(&t,"t.dat",1, 0);
  
  // Integration
  while(n++ < nmax){
    bc1d(f,nx,xoff,1);
    dif1di(f,kx,dt,dx,nx,xoff);
    t+=dt;

    // Output
    if (n % nrec == 0){
      fout( f,"f.dat",nx,1);
      fout(&t,"t.dat",1, 1);
    }
  }
  
  return 0;
}
