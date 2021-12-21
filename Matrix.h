// This file contains the declaration of the Matrix class
#include "Types.h"
#include "Vector.h"
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

  void    setCol(const index j, const Vector x);

  void    show() const;

  void    load(const std::string filename);

  void    write(const std::string filename) const;

  double& operator()(const index i, const index j);

  double  operator()(const index i, const index j) const;

  // Matrix& operator=(const Matrix* A);

};

// here are regular functions for operator overloading

Matrix operator-(const Matrix&);

Matrix operator+(const Matrix&, const Matrix&);

Matrix operator-(const Matrix&, const Matrix&);

Vector operator*(const Matrix&, const Vector&);

std::ostream& operator<<(std::ostream&, const Matrix&);

#endif
