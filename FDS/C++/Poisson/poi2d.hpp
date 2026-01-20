#ifndef _POI2D_HPP_
#define _POI2D_HPP_

#include <iostream>
#include <cmath>
#include "bound.hpp"

double poi2d_ja(double* f, double* rhs,
		double dx, double dy,
		int nx, int ny, int xoff=1, int yoff=1,
		int dnx=0, int dny=0,
		double eps=1e-6, int itmax=128, double omega=2.0/3.0);
double poi2d_gs(double* f, double* rhs,
		double dx, double dy,
		int nx, int ny, int xoff=1, int yoff=1,
		int dnx=0, int dny=0,
		double eps=1e-6, int itmax=128, double omega=1.0);

#endif
