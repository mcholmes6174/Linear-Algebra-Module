#include "Types.h"

#ifndef LSQ
#define LSQ

void decompQR(mat_t& A_OUT, mat_t& Q_OUT);

void solveTriLSQ(const mat_t& Q, const mat_t& R, const vec_t& b, vec_t& x_OUT);

#endif
