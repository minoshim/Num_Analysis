#include "adv1d.hpp"
#include <cmath>

static double median3(double a, double b, double c)
{
  double low=min(b,c);
  double high=max(b,c);
  return max(low,min(a,high));
}

static double pfc_flux(double f0, double a1, double a2, double xi, int sgnv)
{
  const double c0=f0-a2/12.0;
  if (sgnv > 0){
    const double xm=0.5*(1.0-xi);
    const double xl=0.5-xi;
    return (c0+0.5*a1+0.25*a2
	    +4.0*(c0+a1*xm+a2*xm*xm)
	    +c0+a1*xl+a2*xl*xl)/6.0;
  } else{
    const double xm=-0.5*(1.0-xi);
    const double xr=-0.5+xi;
    return (c0-0.5*a1+0.25*a2
	    +4.0*(c0+a1*xm+a2*xm*xm)
	    +c0+a1*xr+a2*xr*xr)/6.0;
  }
}

static void pfc_coef(double* f, int i, double& a1, double& a2)
{
  const double fm=f[i-1];
  const double f0=f[i];
  const double fp=f[i+1];
  const double fmax=max(max(fm,fp),min(2.0*f0-fm,2.0*f0-fp));
  const double fmin=max(0.0,min(min(fm,fp),max(2.0*f0-fm,2.0*f0-fp)));
  const double alpha=1.0/3.0;
  double ap=fp-f0;
  double am=f0-fm;
  const double fpmin=3.0*max(2.0*(f0-fmax),fmin-f0);
  const double fpmax=3.0*min(2.0*(f0-fmin),fmax-f0);
  const double fmmin=-fpmax;
  const double fmmax=-fpmin;

  ap=median3(ap,alpha*fpmin,alpha*fpmax);
  am=median3(am,alpha*fmmin,alpha*fmmax);
  a1=0.5*(ap+am);
  a2=0.5*(ap-am);
}

static void pfc_bounds(double* f, int i, double& fmin, double& fmax)
{
  const double r=2.0/3.0;
  const double fml=f[i-1]+r*(f[i-1]-f[i-2])+(1.0-r)*(f[i]-f[i-1]);
  const double fmr=f[i]+r*(f[i]-f[i+1])+(1.0-r)*(f[i-1]-f[i]);
  const double fpl=f[i]+r*(f[i]-f[i-1])+(1.0-r)*(f[i+1]-f[i]);
  const double fpr=f[i+1]+r*(f[i+1]-f[i+2])+(1.0-r)*(f[i]-f[i+1]);
  const double fmmax=max(max(f[i-1],f[i]),min(fml,fmr));
  const double fmmin=min(min(f[i-1],f[i]),max(fml,fmr));
  const double fpmax=max(max(f[i],f[i+1]),min(fpl,fpr));
  const double fpmin=min(min(f[i],f[i+1]),max(fpl,fpr));

  fmax=max(fmmax,fpmax);
  fmin=max(0.0,min(fmmin,fpmin));
}

static void pfc_correct_coef(double f0, double fmin, double fmax,
			     double a1o, double a2o, double& a1, double& a2)
{
  const double eps=1e-7;
  double ap=a1o+a2o;
  double am=a1o-a2o;
  const double fpmin=3.0*max(2.0*(f0-fmax),fmin-f0);
  const double fpmax=3.0*min(2.0*(f0-fmin),fmax-f0);
  const double fmmin=-fpmax;
  const double fmmax=-fpmin;
  double bp=(ap > 0.0)?ap/(fpmax+eps):ap/(fpmin-eps);
  double bm=(am > 0.0)?am/(fmmax+eps):am/(fmmin-eps);

  bp=min(bp+eps,1.0);
  bm=min(bm+eps,1.0);
  const double apwt=bp/(bp+bm);
  const double amwt=bm/(bp+bm);

  ap=median3(ap,apwt*fpmin,apwt*fpmax);
  am=median3(am,amwt*fmmin,amwt*fmmax);
  a1=0.5*(ap+am);
  a2=0.5*(ap-am);
}

static void wpfc_coef(double* f, int i, double xi, int sgnv, double& a1, double& a2)
{
  const double eps=1e-7;
  const double cwpfc=0.5;
  double fmin,fmax;
  double a1o[3],a2o[3],a1c[3],a2c[3],d[3],gm[3],w[3];

  pfc_bounds(f,i,fmin,fmax);

  a1o[0]=(3.0*f[i]-4.0*f[i-1]+f[i-2])*0.5;
  a2o[0]=(f[i]-2.0*f[i-1]+f[i-2])*0.5;
  a1o[1]=(f[i+1]-f[i-1])*0.5;
  a2o[1]=(f[i+1]-2.0*f[i]+f[i-1])*0.5;
  a1o[2]=(-f[i+2]+4.0*f[i+1]-3.0*f[i])*0.5;
  a2o[2]=(f[i+2]-2.0*f[i+1]+f[i])*0.5;

  for (int k=0;k<3;k++){
    pfc_correct_coef(f[i],fmin,fmax,a1o[k],a2o[k],a1c[k],a2c[k]);
  }

  d[1-sgnv]=(2.0+3.0*xi+xi*xi)/20.0;
  d[1]=(6.0+xi-xi*xi)/10.0;
  d[1+sgnv]=(6.0-5.0*xi+xi*xi)/20.0;

  const double b1=(5.0*a1o[0]+14.0*a1o[1]+5.0*a1o[2])/24.0;
  const double b2=(-a2o[0]+10.0*a2o[1]-a2o[2])/8.0;
  const double b3=(-a1o[0]+2.0*a1o[1]-a1o[2])/6.0;
  const double b4=(a2o[0]-2.0*a2o[1]+a2o[2])/12.0;
  const double dl2=b1*b1/12.0+b2*b2/180.0+b3*b3/448.0+b4*b4/3600.0
    +b1*b3/40.0+b2*b4/420.0;
  double gsum=0.0;

  for (int k=0;k<3;k++){
    const double dl2k=a1c[k]*a1c[k]/12.0+a2c[k]*a2c[k]/180.0;
    gm[k]=d[k]*(cwpfc+sqrt((dl2k+eps)/(dl2+eps)));
    gsum+=gm[k];
  }

  a1=0.0;
  a2=0.0;
  for (int k=0;k<3;k++){
    w[k]=gm[k]/gsum;
    a1+=w[k]*a1c[k];
    a2+=w[k]*a2c[k];
  }
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
    int ic=i-(1+sgnv)/2;
    pfc_coef(f,ic,a1,a2);
    flux[i]=pfc_flux(f[ic],a1,a2,anu,sgnv);
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}

void wpfc(double* f, double v, double dt, double dx, int nx, int xoff)
// Solve 1D advection equation using weighted positive and flux conservative scheme (CFL<1 only)
// f = dependent variable
// v = advection velocity
// dt, dx = time step and grid width
// nx, xoff = number of grid points (including boundary) and boundary grid
// Note: nx-3*xoff corresponds to the number of grid points in computational domain
{
  int i;
  const double nu=v*dt/dx;
  const double anu=(nu >= 0.0)?nu:-nu;
  const int sgnv=(v > 0)?1:-1;
  double flux[nx];

  if (anu == 0.0) return;

  for (i=3;i<nx-2;i++){
    double a1,a2;
    int ic=i-(1+sgnv)/2;
    wpfc_coef(f,ic,anu,sgnv,a1,a2);
    flux[i]=pfc_flux(f[ic],a1,a2,anu,sgnv);
  }
  for (i=xoff;i<nx-xoff;i++) f[i]-=nu*(flux[i+1]-flux[i]);
}
