#include "Constants.h" // for consts::XXX user-defined namespace
#include "Types.h"
#include <string>
#include <vector>

#ifndef LIN_ALG_TOOLS
#define LIN_ALG_TOOLS

namespace tools {

  void catchFileError(const std::string filename, const bool open_failed);

  void catchSingular(const double);

  void checkSymm(const mat_t&);

  vec_t diagPreCond(const mat_t&, const vec_t&);

  vec_t getError(const mat_t&, const vec_t&, const vec_t&);

  double getNorm(const vec_t&);

  double innerProd(const vec_t&, const vec_t&);

  mat_t matMul(const mat_t&, const mat_t&);

  vec_t makeVecCopy(const vec_t&);

  vec_t matVecMul(const mat_t&, const vec_t&);

  mat_t readMatrix(const std::string);

  int readNthVal(const std::string, const int);

  vec_t readVector(const std::string);

  void showMatrix(const mat_t&);

  void showVector(const vec_t&);

  void swapRows(mat_t&, vec_t&, const int, const int);

  void writeMatrix(const std::string, const mat_t&);

  void writeVector(const std::string, const vec_t&);

}

#endif
