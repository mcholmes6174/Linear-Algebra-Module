#include "Types.h"

#ifndef GMRES
#define GMRES

void applyGMRES(mat_t& A_OUT, const vec_t& b, vec_t& x_OUT, const size_t n);

#endif
