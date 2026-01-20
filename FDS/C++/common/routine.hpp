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
#include "bound.hpp"

template <typename T>
int fout(const T* x, const char s[100], int nx, bool flag)
// Output ascii file
{
  FILE* fil;
  const char* fmt=nullptr;

  if constexpr (std::is_same_v<T, double>) {
    fmt = "%.15f\n";
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

#endif
