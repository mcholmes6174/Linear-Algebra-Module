// This file contains the declaration of the Matrix class.
#include "Types.h"

#ifndef MATRIX
#define MATRIX

class Matrix {

private:
  index m_row{};
  index m_col{};
  mat_t m_mat{};

public:
  Matrix(const index m=1, const index n=1)
    : m_row{m}, m_col{n}, m_mat(m,vec_t(n)) {
    // default constructor for our Matrix class
  }

  void   setSize(const index m, const index n);

  index  getSize(const int dim);

  void   setVal(const index i, const index j, const double a_ij);

  double getVal(const index i, const index j);

  void   show();
};

#endif
