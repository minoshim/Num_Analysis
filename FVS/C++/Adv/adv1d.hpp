#ifndef _ADV1D_HPP_
#define _ADV1D_HPP_

void fv3rd(double* f, double v, double dt, double dx, int nx, int xoff=2);
void muscl(double* f, double v, double dt, double dx, int nx, int xoff=2);

inline double max(double a, double b)
{
  return (a >= b)?a:b;
}
inline double min(double a, double b)
{
  return (a >= b)?b:a;
}
inline double minmod(double a, double b)
{
  return min(max(a,b),0)+max(min(a,b),0);
}
inline double muscl_mm_cal_f(double fu, double f0, double fd)
{
  // Minmod
  return f0-0.5*minmod(fu-f0,f0-fd);
}
inline double muscl_mc_cal_f(double fu, double f0, double fd)
{
  // Monotonized central
  double alpha=2.0;
  return f0-0.5*minmod(alpha*minmod(fu-f0,f0-fd),0.5*(fu-fd));
}
inline double muscl_vl_cal_f(double fu, double f0, double fd)
{
  // van Leer
  double d[3]={fu-f0,f0-fd};
  d[2]=d[0]*d[1];
  double df=(d[2]>0)?(d[2]/(d[0]+d[1])):0;
  return f0-df;
}

#endif
