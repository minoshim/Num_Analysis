#ifndef _DIF1D_HPP_
#define _DIF1D_HPP_

#include <cmath>

void dif1de(double* f, double kx, double dt, double dx, int nx, int xoff=1);
void dif1di(double* f, double kx, double dt, double dx, int nx, int xoff=1,
	    double alpha=1.0);
void trdiag(double* s,
	    const double* a, const double* b, const double* c, const double* d,
	    int nx, int xoff=1);

#endif
