#ifndef _ROUTINE_HPP_
#define _ROUTINE_HPP_

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>
#include <cstdlib>
#include <cmath>
#include <time.h>
#include <random>
#include "bound.hpp"

template <typename T>
int fout(const T* x, const char s[100], int nx, bool flag)
// Output ascii file
{
  FILE* fil;
  const char* fmt=nullptr;

  if constexpr (std::is_same_v<T, double>) {
    fmt = "%.9f\n";
  } else if constexpr (std::is_same_v<T, int>) {
    fmt = "%d\n";
  }

  if ( (fil=fopen(s,(flag)?"a":"w")) != NULL){
    for (int i=0;i<nx;i++){
      fprintf(fil,fmt,x[i]);
    }
    fclose(fil);
  } else{
    perror("Error to open file");
    return -1;
  }
  return 0;
}

inline double rand_noise(const double *params, unsigned seed)
{
  // Return uniform random distribution using Mersenne Twister
  // params[0] = mean
  // params[1] = width, giving params[0]+-params[1]=max/min
  static int r_flag=0;
  static std::mt19937_64 engine;
  if (r_flag == 0){
    engine.seed(seed);
    r_flag=1;
  }
  std::uniform_real_distribution<> dist(params[0]-params[1],params[0]+params[1]);
  return dist(engine);
}

#endif
