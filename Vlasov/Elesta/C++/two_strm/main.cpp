#include "cnst.hpp"
#include "routine.hpp"
#include "adv1d.hpp"

int main(void)
{
  int i,j,n;
  const int xoff=XOFF;
  const int voff=VOFF;
  const int offs[]={xoff,voff};
  const int nx=XMESH+2*xoff;
  const int nv=VMESH+2*voff;
  const int nd=nx*nv;
  const int nrec=NREC;
  const int nmax=NMAX;
  const double lx=4.0*M_PI;		// Domain size in X, [-lx/2,lx/2]
  const double lv=4.5;		// Half domain size in V, [-lv,lv]
  const double dx=lx/XMESH;	// Grid spacing in X
  const double dv=2.0*lv/VMESH; 	// Grid spacing in V
  const double dt=fabs(CFL*dx/lv); // Time step
  const double q=-1.0;		   // Normalized charge for electrons 
  double t=0.0;
  double x[nx],v[nv];
  double* g=new double[nx]();
  double* drho=new double[nx]();
  double* f=new double[nd]();

  // Simulation paramters
  const double amp=0.05;	// Perturbation amplitude
  const double kk=0.5;	// Perturbation wavenumber

  // Initialize
  for (i=0;i<nx;i++){
    x[i]=(i+0.5-xoff)*dx-0.5*lx;
  }
  for (j=0;j<nv;j++){
    v[j]=(j+0.5-voff)*dv-lv;
  }
  for (j=0;j<nv;j++){
    double ftmp=(v[j]*v[j])*gaussian(v[j],0.0,1.0);
    for (i=0;i<nx;i++){
      f[nx*j+i]=ftmp*(1.0+amp*cos(kk*x[i]));
    }
  }

  // Output
  fout(offs,"offs.dat",2,0);
  fout(x,"x.dat",nx,0);
  fout(v,"v.dat",nv,0);
  fout(f,"f.dat",nd,0);
  fout(g,"g.dat",nx,0);
  fout(&t,"t.dat",1, 0);
  
  // Time integration
  while(n++ < nmax){

    // Advection in X (half)
    bc2d(f,nx,nv,xoff,voff,0,999); // Periodic in x, nothing to do in v
    for (j=voff;j<nv-voff;j++){
      csl3rd(&f[nx*j],v[j],0.5*dt,dx,nx,xoff);
    }
    
    // Calculate field
    {
      // Density fluctuation (integration and zero-mean)
      double dmean=0.0;
      for (i=xoff;i<nx-xoff;i++){
	drho[i]=0.0;
	for (j=voff;j<nv-voff;j++){
	  drho[i]+=f[nx*j+i]*dv;
	}
	dmean+=drho[i]*dx;
      }
      dmean/=lx;
      for (i=xoff;i<nx-xoff;i++){
	drho[i]-=dmean;
      }

      // Electric field (integration and zero-mean)
      g[xoff]=0.0;
      for (i=xoff+1;i<nx-xoff+1;i++){
	g[i]=g[i-1]+q*drho[i-1]*dx;
      }
      double gmean=0.0;
      for (i=xoff;i<nx-xoff;i++){
	gmean+=0.5*(g[i]+g[i+1])*dx;
      }
      gmean/=lx;
      for (i=xoff;i<nx-xoff+1;i++){
	g[i]-=gmean;
      }
    }
    
    // Advection in V (full)
    for (i=xoff;i<nx-xoff;i++){
      double ftmp[nv];
      for (j=0   ;j<nv     ;j++) ftmp[j]=f[nx*j+i];
      csl3rd(ftmp,q*0.5*(g[i]+g[i+1]),dt,dv,nv,voff);
      for (j=voff;j<nv-voff;j++) f[nx*j+i]=ftmp[j];
    }
    
    // Advection in X (half)
    bc2d(f,nx,nv,xoff,voff,0,999); // Periodic in x, nothing to do in v
    for (j=voff;j<nv-voff;j++){
      csl3rd(&f[nx*j],v[j],0.5*dt,dx,nx,xoff);
    }

    t+=dt;
    
    // Output
    if (n % nrec == 0){
      fout(f,"f.dat",nd,1);
      fout(g,"g.dat",nx,1);
      fout(&t,"t.dat",1, 1);
    }
  }

  
  delete[] g;
  delete[] drho;
  delete[] f;
  return 0;
}
