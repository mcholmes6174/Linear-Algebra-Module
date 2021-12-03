#include "Types.h"

#ifndef CONJUGATE_GRAD
#define CONJUGATE_GRAD

void basicCG(       const mat_t& A, const vec_t& b, vec_t& x);

void smartCG(       const mat_t& A, const vec_t& b, vec_t& x);

void smartPreCondCG(const mat_t& A, const vec_t& b, vec_t& x);

#endif
