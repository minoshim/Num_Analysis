#include "adv1d.hpp"

static double median3(double a, double b, double c)
{
  double low=min(b,c);
  double high=max(b,c);
  return max(low,min(a,high));
}

static double pfc_poly(double f0, double a1, double a2, double x)
{
  return f0-a2/12.0+a1*x+a2*x*x;
}

static double pfc_flux(double f0, double a1, double a2, double xi, int sgnv)
{
  if (sgnv > 0){
    return xi*(pfc_poly(f0,a1,a2,0.5)
	       +4.0*pfc_poly(f0,a1,a2,0.5*(1.0-xi))
	       +pfc_poly(f0,a1,a2,0.5-xi))/6.0;
  } else{
    return xi*(pfc_poly(f0,a1,a2,-0.5)
	       +4.0*pfc_poly(f0,a1,a2,-0.5*(1.0-xi))
	       +pfc_poly(f0,a1,a2,-0.5+xi))/6.0;
  }
}

static void pfc_coef(double* f, int i, double& a1, double& a2)
{
  const double fmax=max(max(f[i-1],f[i+1]),min(2.0*f[i]-f[i-1],2.0*f[i]-f[i+1]));
  const double fmin=max(0.0,min(min(f[i-1],f[i+1]),max(2.0*f[i]-f[i-1],2.0*f[i]-f[i+1])));
  const double alpha=1.0/3.0;
  double ap=f[i+1]-f[i];
  double am=f[i]-f[i-1];
  const double fpmin=3.0*max(2.0*(f[i]-fmax),fmin-f[i]);
  const double fpmax=3.0*min(2.0*(f[i]-fmin),fmax-f[i]);
  const double fmmin=3.0*max(2.0*(fmin-f[i]),f[i]-fmax);
  const double fmmax=3.0*min(2.0*(fmax-f[i]),f[i]-fmin);

  ap=median3(ap,alpha*fpmin,alpha*fpmax);
  am=median3(am,alpha*fmmin,alpha*fmmax);
  a1=0.5*(ap+am);
  a2=0.5*(ap-am);
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

void muscl(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using MUSCL scheme
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
  double (*msl_func)(double,double,double)=muscl_mc_cal_f;
    
  for (i=2;i<nx-1;i++) flux[i]=0.5*(+(1+sgnv)*msl_func(f[i-2],f[i-1],f[i  ])
				    +(1-sgnv)*msl_func(f[i+1],f[i  ],f[i-1]));
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}

void csl3rd(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using conservative semi-Lagrange scheme (CFL<1 only)
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

  for (i=2;i<nx-1;i++){
    double c0[2]={(-f[i-2]+5*f[i-1]+2*f[i  ])/6.0,
		  (-f[i+1]+5*f[i  ]+2*f[i-1])/6.0};
    double c1[2]={f[i]-f[i-1],f[i]-f[i-1]};
    double c2[2]={(f[i-2]-2*f[i-1]+f[i  ])*0.5,
		  (f[i+1]-2*f[i  ]+f[i-1])*0.5};
    double ft[2]={(c0[0]+nu*(-c1[0]*0.5+c2[0]*nu/3.0)),
		  (c0[1]+nu*(-c1[1]*0.5+c2[1]*nu/3.0))};
    flux[i]=0.5*(+(1+sgnv)*ft[0]
		 +(1-sgnv)*ft[1]);
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}

void cslmsl(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using conservative semi-Lagrange-MUSCL scheme (CFL<1 only)
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

  for (i=2;i<nx-1;i++){
    double c0[2]={(-f[i-2]+5*f[i-1]+2*f[i  ])/6.0,
		  (-f[i+1]+5*f[i  ]+2*f[i-1])/6.0};
    double c1[2]={f[i]-f[i-1],f[i]-f[i-1]};
    double c2[2]={(f[i-2]-2*f[i-1]+f[i  ])*0.5,
		  (f[i+1]-2*f[i  ]+f[i-1])*0.5};
    double ft[2]={(c0[0]+nu*(-c1[0]*0.5+c2[0]*nu/3.0)),
		  (c0[1]+nu*(-c1[1]*0.5+c2[1]*nu/3.0))};
    flux[i]=0.5*(+(1+sgnv)*(f[i-1]+minmod(f[i-1]-f[i-2],f[i  ]-f[i-1],ft[0]-f[i-1]))
		 +(1-sgnv)*(f[i  ]-minmod(f[i  ]-f[i-1],f[i+1]-f[i  ],f[i  ]-ft[1])));
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}

void pfc(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using positive and flux conservative scheme (CFL<1 only)
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-2*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  const double anu=(nu >= 0.0)?nu:-nu;
  const int sgnv=(v > 0)?1:-1;
  double flux[nx];

  if (anu == 0.0) return;

  for (i=2;i<nx-1;i++){
    double a1,a2;
    int ic=(sgnv > 0)?i-1:i;
    pfc_coef(f,ic,a1,a2);
    flux[i]=pfc_flux(f[ic],a1,a2,anu,sgnv)/anu;
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}
