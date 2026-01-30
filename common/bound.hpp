#ifndef _BOUND_HPP_
#define _BOUND_HPP_

#include <cmath>

void bc1d(double *f, int nx, int xoff, int dnx=0);
void bc2d(double *f, int nx, int ny, int xoff, int yoff, int dnx=0, int dny=0);

#endif
