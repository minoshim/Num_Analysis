#ifndef _BGS1D_HPP_
#define _BGS1D_HPP_

void bgs_fd(double* f, double dt, double dx, int nx, int xoff);
void bgs_fv(double* f, double dt, double dx, int nx, int xoff);

#endif
