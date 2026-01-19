#include "poi2d.hpp"

double poi2d_ja(double* f, double* rhs,
		double dx, double dy,
		int nx, int ny, int xoff, int yoff,
		int dnx, int dny,
		double eps, int itmax, double omega)
{
  // Solve 2D poisson equation d2f/dx2+d2f/dy2 = -rhs using weighted Jacobi method
  // dx and dy are constant grid width
  // omega(<=1) is optimization paramter, typically 2/3

  int i,j,ss;
  int cnt=0;
  const double eps2=eps*eps;
  const double dx2=dx*dx,dy2=dy*dy,dxdy2=dx2*dy2;
  const double denom=1./(2.0*(dx2+dy2));
  double anorm,anormf=0.0;
  double* fold=new double[nx*ny];

  // Boundary condition

  for (j=yoff;j<ny-yoff;j++){
    for (i=xoff;i<nx-xoff;i++){
      double resid=dxdy2*rhs[nx*j+i]*denom;
      anormf+=resid*resid;
    }
  }
  anormf=(anormf <= eps2)?eps2:anormf;
  
  do{
    anorm=0.0;

    for (ss=0;ss<nx*ny;ss++) fold[ss]=f[ss];
    
    for (j=yoff;j<ny-yoff;j++){
      for (i=xoff;i<nx-xoff;i++){
	ss=nx*j+i;
	double resid=(+dy2*(fold[ss- 1]+fold[ss+ 1])
		      +dx2*(fold[ss-nx]+fold[ss+nx])
		      +dxdy2*rhs[ss])*denom-fold[ss];
	resid*=omega;
	anorm+=resid*resid;
	f[ss]+=resid;
      }
    }
    // Boundary condition
    
  } while( (anorm > eps2*anormf) && (++cnt < itmax) );

  delete[] fold;
  return sqrt(anorm/anormf);
}
  
double poi2d_gs(double* f, double* rhs,
		double dx, double dy,
		int nx, int ny, int xoff, int yoff,
		int dnx, int dny,
		double eps, int itmax, double omega)
  // Solve 2D poisson equation d2f/dx2+d2f/dy2 = -rhs using red-black Gauss-Seidel method
  // dx and dy are constant grid width
  // omega(<=2) is optimization paramter for SOR method

  int i,j,ss;
  int cnt=0;
  const double eps2=eps*eps;
  const double dx2=dx*dx,dy2=dy*dy,dxdy2=dx2*dy2;
  const double denom=1./(2.0*(dx2+dy2));
  double anorm,anormf=0.0;

  // Boundary condition

  for (j=yoff;j<ny-yoff;j++){
    for (i=xoff;i<nx-xoff;i++){
      double resid=dxdy2*rhs[nx*j+i]*denom;
      anormf+=resid*resid;
    }
  }
  anormf=(anormf <= eps2)?eps2:anormf;
  
  do{
    anorm=0.0;
    for (int pass=0;pass<2;pass++){	// Red-Black
      for (j=yoff;j<ny-yoff;j++){
	int isw=xoff+1-((j+pass) % 2);
	for (i=isw;i<nx-xoff;i+=2){
	  ss=nx*j+i;
	  double resid=(+dy2*(f[ss- 1]+f[ss+ 1])
			+dx2*(f[ss-nx]+f[ss+nx])
			+dxdy2*rhs[ss])*denom-f[ss];
	  resid*=omega;
	  anorm+=resid*resid;
	  f[ss]+=resid;
	}
      }
    }
    // Boundary condition
    
  } while( (anorm > eps2*anormf) && (++cnt < itmax) );

  return sqrt(anorm/anormf);
{
}

