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
#include "ftcs.hpp"

int fout(const double* x, const char s[100], int nx, bool flag)
// Output ascii file
{
  FILE* fil;
  if ( (fil=fopen(s,(flag)?"a":"w")) != NULL){
    for (int i=0;i<nx;i++){
      fprintf(fil,"%.15f\n",x[i]);
    }
    fclose(fil);
  } else{
    perror("Error to open file");
    return -1;
  }
  return 0;
}

#endif
