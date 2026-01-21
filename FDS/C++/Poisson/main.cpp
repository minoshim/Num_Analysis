#include "cnst.hpp"
#include "routine.hpp"
#include "poi2d.hpp"

int main(void)
{
  int i,j,n;
  const int xoff=XOFF;
  const int yoff=YOFF;
  const int offs[]={xoff,yoff};
  const int nx=XMESH+2*xoff;
  const int ny=YMESH+2*yoff;
  const int nd=nx*ny;
  const double lx=2*M_PI;	// Domain size
  const double ly=2*M_PI;	// Domain size
  const double dx=lx/XMESH;	// Grid spacing
  const double dy=ly/YMESH;	// Grid spacing
  const double params[2]={0.0,1.0}; // Parameter for uniform random dist.
  const unsigned seed=12345;
  double *x,*y,*f,*rhs;
  x=new double[nx]();
  y=new double[ny]();
  f=new double[nd]();
  rhs=new double[nd]();

  // Initialize
  for (i=0;i<nx;i++){
    x[i]=(i+0.5-xoff)*dx-0.5*lx;
  }
  for (j=0;j<ny;j++){
    y[j]=(j+0.5-yoff)*dy-0.5*ly;
  }
  for (j=0;j<ny;j++){
    for (i=0;i<nx;i++){
      int ss=nx*j+i;
      f[ss]=0.0;
      // rhs[ss]=rand_noise(params,seed);
      rhs[ss]=sin(2*x[i])*sin(3*y[j]);
    }
  }
  bc2d(rhs,nx,ny,xoff,yoff,0,0);
  
  // Solve Poisson equation
  poi2d_ja(f,rhs,dx,dy,nx,ny,xoff,yoff,0,0,1e-6,8192);
  
  // Output
  fout(offs,"offs.dat",2,0);
  fout(x,"x.dat",nx,0);
  fout(y,"y.dat",ny,0);
  fout(f,"f.dat",nd,0);
  fout(rhs,"rhs.dat",nd,0);
  
  delete[] x;
  delete[] y;
  delete[] f;
  delete[] rhs;
}
