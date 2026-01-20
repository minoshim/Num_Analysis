#include "bound.hpp"

void bc1d(double *f, int nx, int xoff, int dnx)
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

void bc2d(double *f, int nx, int ny, int xoff, int yoff, int dnx, int dny)
// 2D boundary condition
// dn: 0 for periodic, -1 for fix, +1 for free
{
  int i,j;
  if (dnx == 0){
    // Periodic
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f[nx*j+(nx-1-i)]=f[nx*j+(2*xoff-1-i)];
	f[nx*j+i]=f[nx*j+(nx-2*xoff+i)];
      }
    }
  } else if (abs(dnx) == 1){
    // Fix or free
    for (j=0;j<ny;j++){
      for (i=0;i<xoff;i++){
	f[nx*j+i]=dnx*f[nx*j+(2*xoff-1)-i];
      }
      for (i=0;i<xoff;i++){
	f[nx*j+(nx-1-i)]=dnx*f[nx*j+(nx-2*xoff)+i];
      }
    }
  } else{
    // Nothing to do
  }
  
  if (dny == 0){
    // Periodic
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*(ny-1-j)+i]=f[nx*(2*yoff-1-j)+i];
	f[nx*j+i]=f[nx*(ny-2*yoff+j)+i];
      }
    }
  } else if (abs(dny) == 1){
    // Fix or free
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*j+i]=dny*f[nx*(2*yoff-1-j)+i];	
      }
    }
    for (j=0;j<yoff;j++){
      for (i=0;i<nx;i++){
	f[nx*(ny-1-j)+i]=dny*f[nx*(ny-2*yoff+j)+i];
      }
    }
  } else {
    // Nothing to do
  }
}
