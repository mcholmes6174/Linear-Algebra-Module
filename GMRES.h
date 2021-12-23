#include "Matrix.h"
#include "Vector.h"

#ifndef GMRES
#define GMRES

void applyGMRES(const Matrix& A_OUT, const Vector& b_OUT, Vector& x_OUT);

#endif
