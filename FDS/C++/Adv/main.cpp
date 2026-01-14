#include "cnst.hpp"
#include "routine.hpp"

int main(void)
{
  int i,n;
  const int nx=XMESH+2*XOFF;
  const int nrec=10;
  const int nmax=nrec*100;
  const double lx=1.0;
  const double dx=lx/XMESH;
  const double v=1.0;
  const double dt=fabs(CFL*dx/v);
  double t=0.0;
  double *x,*f;
  x=new double[nx]();
  f=new double[nx]();

  // Initialize
  for (i=0;i<nx;i++){
    x[i]=(i+0.5-XOFF)*dx;
    f[i]=(i > nx/4 && i < 3*nx/4)?1.0:0.0;
  }

  // Output
  fout( x,"x.dat",nx,0);
  fout( f,"f.dat",nx,0);
  fout(&t,"t.dat",1, 0);

  // Integration
  while(n++ < nmax){
    bound(f,nx,XOFF,0);
    ftcs(f,v,dt,dx,nx,XOFF);
    t+=dt;

    // Output
    if (n % nrec == 0){
      fout( f,"f.dat",nx,1);
      fout(&t,"t.dat",1, 1);
    }
  }
  
  return 0;
}
