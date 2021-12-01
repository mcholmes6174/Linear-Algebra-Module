#include "Constants.h"
#include "Types.h"
#include <vector>

#ifndef GAUSSIAN_ELIM
#define GAUSSIAN_ELIM

void pivotElim(mat_t& A, vec_t& b);

void backSub(const mat_t& U, const vec_t& y, vec_t& x);

#endif
