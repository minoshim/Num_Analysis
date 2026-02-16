#include "cnst.hpp"
#include "routine.hpp"
#include "adv1d.hpp"

int main(void)
{
  int i,j,n;
  const int xoff=XOFF;
  const int voff=VOFF;
  const int nx=XMESH+2*xoff;
  const int nv=VMESH+2*voff;
  const int nrec=NREC;
  const int nmax=NMAX;
  const double lx=1.0;		// Domain size in X
  const double lv=2.0;		// Domain size in V
  const double dx=lx/XMESH;	// Grid spacing in X
  const double dv=lv/VMESH; 	// Grid spacing in V
  const double dt=fabs(CFL*dx/(0.5*lv)); // Time step
  double t=0.0;
  double x[nx],v[nv];
  double* f=new double[nx*nv]();

  delete[] f;
  return 0;
}
