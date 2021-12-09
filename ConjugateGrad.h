#include "Matrix.h"
#include "Vector.h"

#ifndef CONJUGATE_GRAD
#define CONJUGATE_GRAD

void basicCG(       const Matrix& A, const Vector& b, Vector& x);

void smartCG(       const Matrix& A, const Vector& b, Vector& x);

void smartPreCondCG(const Matrix& A, const Vector& b, Vector& x);

#endif
