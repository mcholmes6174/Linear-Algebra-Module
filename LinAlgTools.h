#include "Constants.h" // for consts::XXX user-defined namespace
#include <string>

#ifndef LIN_ALG_TOOLS
#define LIN_ALG_TOOLS

namespace tools {

  void catchSingular(const double value);

  void checkSymm(const double A[][consts::n], const int m=consts::m);

  void diagPreCond(const double M[][consts::n], const double r[], double z[],
                   const int m=consts::m);

  void getError(const double A[][consts::n], const double   x[],
                const double b[],                  double err[],
                const int m=consts::m, const int n=consts::n);

  double getNorm(const double x[], const int m=consts::m);

  double innerProd(const double x[], const double y[], const int m=consts::m);

  void matMul(const double A[][consts::n], const double B[][consts::n],
                    double C[][consts::n], const int m=consts::m,
                                           const int n=consts::n);

  void makeVecCopy(double x_copy[], const double x[], const int m=consts::m);

  void matVecMul(const double A[][consts::n], const double x[], double y[],
                 const int m=consts::m, const int n=consts::n);

  void readMatrix(const std::string filename, double A[][consts::n],
                  const int m=consts::m, const int n=consts::n);

  void showMatrix(const double A[][consts::n], const int m=consts::m,
                                               const int n=consts::n);

  void showVector(const double x[], const int m=consts::m);

  void swapRows(double A[][consts::n], double b[],
                const int row1, const int row2, const int n=consts::n);

  void writeVector(const std::string filename, const double x[],
                                               const int m=consts::m);

}

#endif
