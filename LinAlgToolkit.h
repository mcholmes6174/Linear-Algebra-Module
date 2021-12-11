#include "Constants.h"
#include "Matrix.h"
#include "Types.h"
#include "Vector.h"
#include <string>
#include <vector>

#ifndef LIN_ALG_TLK
#define LIN_ALG_TLK

/*
 * References used as function parameters (e.g., mat_t&) allow us to avoid
 * making a copy of the argument being passed, and allows the argument to be
 * modified within the function. If a reference is passed as a constant (e.g.,
 * const mat_t&), then the argument is still accessed directly, but is
 * guaranteed not to be modified within the function.
 *
 * Functions that would otherwise return vec_t or mat_t types utilize non-const
 * reference parameters for performance reasons. These parameters are labeled
 * with the suffix _OUT. All other functions return by value if not void.
*/

namespace tlk {

  void   catchFileError(const std::string, const bool);

  void   catchSingular(const double);

  void   checkSymm(const Matrix&);

  void   diagPreCond(const Matrix&, const Vector&, Vector& z_OUT);

  void   getError(const Matrix&, const Vector&, const Vector&, Vector& err_OUT);

  double innerProd(const Vector&, const Vector&);

  void   matMul(const Matrix&, const Matrix&, Matrix& C_OUT);

  // we return by value here b/c we want a copy
  Vector makeVecCopy(const Vector&, const double scalar=1.0);

  void   matVecMul(const Matrix&, const Vector&, Vector& y_OUT);

  void   getCol(const index j, Matrix& A, Vector& x);

  index  readNthVal(const std::string, const int);

  void   swapRows(Matrix& A_OUT, Vector& b_OUT, const index, const index);

  // we return by value here to avoid resizing and additional transpositions
  Matrix transpose(const Matrix&);

  // we return by value here to increase readability
  Vector updateVector(const Vector&, const double, const Vector&);

}

#endif
