// This file contains the declaration of the Matrix class.
#include "Types.h"
#include <string>

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

  void    resize(const index m, const index n);

  index   size(const int dim) const;

  double& set(const index i, const index j);

  double  get(const index i, const index j) const;

  void    show() const;

  void    load(const std::string filename);

  void    write(const std::string filename) const;

};

#endif
