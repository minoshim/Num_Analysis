#include "bound.hpp"

void bound(double *f, int nx, int xoff, int dnx)
// 1D boundary condition
// dn: 0 for periodic, -1 for fix, +1 for free
{
  int i;
  if (dnx == 0){
    // Periodic
    for (i=0;i<xoff;i++){
      f[i]=f[nx-2*xoff+i];      
      f[nx-1-i]=f[2*xoff-1-i];
    }
  } else if (abs(dnx) == 1){
    // Fix or free
    for (i=0;i<xoff;i++){
      f[i]=dnx*f[2*xoff-1-i];
      f[nx-1-i]=dnx*f[nx-2*xoff+i];
    }
  } else{
    // Nothing to do
  }
}

