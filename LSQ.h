#include "Types.h"

#ifndef LSQ
#define LSQ

void decompQR(mat_t& A_OUT, vec_t& b_OUT, vec_t& Q_OUT);

void solveTriLSQ(const mat_t& A_OUT, const vec_t& b, vec_t& x_OUT);

#endif
