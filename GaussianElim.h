#include "Matrix.h"
#include "Vector.h"

#ifndef GAUSSIAN_ELIM
#define GAUSSIAN_ELIM

void pivotElim(Matrix& A, Vector& b);

void backSub(const Matrix& U, const Vector& y, Vector& x);

#endif
