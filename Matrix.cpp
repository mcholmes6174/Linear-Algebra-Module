// This file contains the implementations of the member functions for the
// Matrix class
#include "Constants.h"
#include "Matrix.h"
#include "Types.h"
#include <cassert>
#include <iomanip>
#include <iostream>

void Matrix::setSize(const index m, const index n) {
  // to reset the size of the matrix
  m_row = m;
  m_col = n;
  m_mat.resize(m_row,vec_t(m_col));
}

index Matrix::getSize(const int dim) {
  // to get the size of the matrix
  if      (dim == 0) return m_mat.size();
  else if (dim == 1) return m_mat[0].size();
  else std::cerr << "Must choose dim = 0 or dim = 1 in Matrix.getSize().";
  assert(dim==0 || dim==1);
}

void Matrix::setVal(const index i, const index j, const double a_ij) {
  // to set the value of an individual entry
  assert(i < m_row && j < m_col);
  m_mat[i][j] = a_ij;
}

double Matrix::getVal(const index i, const index j) {
  // to get the value stored at index i,j
  assert(i < m_row && j < m_col);
  return m_mat[i][j];
}

void Matrix::show() {
  // to print the matrix to the screen
  using std::cout;
  cout << std::setprecision(consts::dispPrec);
  for (index i{}; i < m_row; ++i) {
    cout << '\n';
    for (auto a_ij : m_mat[i]) {
      cout << std::setw(consts::dispPrec+2)
           << std::right << std::showpoint << a_ij << '\t';
    }
  }
  cout << '\n';
}
