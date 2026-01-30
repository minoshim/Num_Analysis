#ifndef _ADV1D_HPP_
#define _ADV1D_HPP_

void ftcs(double* f, double v, double dt, double dx, int nx, int xoff=1);
void upwd(double* f, double v, double dt, double dx, int nx, int xoff=1);
void lawe(double* f, double v, double dt, double dx, int nx, int xoff=1);
void fv3rd(double* f, double v, double dt, double dx, int nx, int xoff=2);

#endif
