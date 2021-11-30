#include "Constants.h"

#ifndef GAUSSIAN_ELIM
#define GAUSSIAN_ELIM

void pivotElim(double mat_A[][consts::n], double vec_b[],
               const int m=consts::m, const int n=consts::n);

void backSub(const double mat_U[][consts::n],
             const double vec_y[], double vec_x[],
             const int m=consts::m, const int n=consts::n);

#endif
