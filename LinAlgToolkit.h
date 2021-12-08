#include "Classes.h"
#include "Constants.h"
#include "Types.h"
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

  void   checkSymm(const mat_t&);

  void   diagPreCond(const mat_t&, const vec_t&, vec_t& z_OUT);

  void   getError(const mat_t&, const vec_t&, const vec_t&, vec_t& err_OUT);

  double getNorm(const vec_t&);

  double innerProd(const vec_t&, const vec_t&);

  double innerProd2(Vector&, Vector&);

  void   matMul(const mat_t&, const mat_t&, mat_t& C_OUT);

  // we return by value here b/c we want a copy
  vec_t  makeVecCopy(const vec_t&, const double scalar=1.0);

  void   matVecMul(const mat_t&, const vec_t&, vec_t& y_OUT);

  void   normalizeVec(vec_t& x_OUT);

  void   nthColumn(const std::string, mat_t&, vec_t&, const index);

  void   readMatrix(const std::string, mat_t& A_OUT);

  index  readNthVal(const std::string, const int);

  void   readVector(const std::string, vec_t& x_OUT);

  void   showMatrix(const mat_t&);

  void   showVector(const vec_t&);

  void   swapRows(mat_t& A_OUT, vec_t& b_OUT, const index, const index);

  // we return by value here to avoid resizing and additional transpositions
  mat_t  transpose(const mat_t&);

  // we return by value here to increase readability
  vec_t  updateVector(const vec_t&, const double, const vec_t&);

  void   writeMatrix(const std::string, const mat_t&);

  void   writeVector(const std::string, const vec_t&);

}

#endif
