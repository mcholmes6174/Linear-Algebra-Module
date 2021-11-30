#include "Constants.h"

#ifndef CONJUGATE_GRAD
#define CONJUGATE_GRAD

void basicCG(       const double A[][consts::n], const double b[], double x[]);

void smartCG(       const double A[][consts::n], const double b[], double x[]);

void smartPreCondCG(const double A[][consts::n], const double b[], double x[]);

#endif
