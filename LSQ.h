#include "Matrix.h"
#include "Vector.h"

#ifndef LSQ
#define LSQ

void decompQR(Matrix& A_OUT, Vector& b_OUT, Vector& Q_OUT);

void solveTriLSQ(const Matrix& A_OUT, const Vector& b, Vector& x_OUT);

#endif
