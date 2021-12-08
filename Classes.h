// This file contains all the user-defined classes that are used throughout the
// program.
#include "Constants.h"
#include "Types.h"
#include <cmath>
#include <iomanip>
#include <iostream>

#ifndef CLASSES
#define CLASSES

class Vector {
// private members can only be accessed by other members of Vector
private:
  index m_row{};  // a std::size_t datatype
  vec_t m_vec{};  // a std::vector datatype

// public members can be accessed by anyone
public:
  Vector() = default; // a blank default constructor

  Vector(index m = 1) { // a non-default constructor
    m_row = m;
    m_vec.resize(m_row);
  }

  void setSize(const index m) {
    // to set the length of the vector
    m_row = m;
    m_vec.resize(m_row);
  }

  index getSize() {
    // to get the length of the vector
    return m_vec.size();
  }

  void setVal(const index i, const double v_i) {
    // to set the value of a particular entry
    if (i >= m_row) return;
    m_vec[i] = v_i;
  }

  double getVal(const index i) {
    // to get the value stored at index i
    if (i >= m_row) {
      std::cout << "Error: Index out of range in Vector.getVal(). Returning 0.";
      return 0;
    }
    return m_vec[i];
  }

  void show() {
    // to print the vector to the screen
    using std::cout;
    cout << std::setprecision(consts::dispPrec);
    for (auto v_i : m_vec) {
      cout << '\n' << std::setw(consts::dispPrec+2)
           << std::right << std::showpoint << v_i;
    }
    cout << '\n';
  }

  double getNorm(const double p = 2.0) {
    // to return the p-norm (2 by default) norm of the vector
    double norm{};
    for (auto v_i : m_vec) {
      norm += std::pow(v_i, p);
    }
    return    std::pow(norm, 1.0/p);
  }
};

//****************************************************************************//

class Matrix {

private:
  index m_row{};
  index m_col{};
  mat_t m_mat{};

public:
  void setSize(const index m, const index n) {
    // to set the size of the matrix
    m_row = m;
    m_col = n;
    m_mat.resize(m_row,vec_t(m_col));
  }

  index getSize(const int dim) {
    // to get the size of the matrix
    if      (dim == 0) return m_mat.size();
    else if (dim == 1) return m_mat[0].size();
    else std::cout << "Must choose dim = 0 or dim = 1 in Matrix.getSize().";
  }

  void setVal(const index i, const index j, const double a_ij) {
    // to set the value of an individual entry
    if (i >= m_row || j >= m_col) return;
    m_mat[i][j] = a_ij;
  }

  double getVal(const index i, const index j) {
    // to get the value stored at index i,j
    if (i >= m_row) {
      std::cout << "Error: Index i out of range in Matirx.getVal(). Returning 0.";
      return 0;
    }
    else if (j >= m_col){
      std::cout << "Error: Index j out of range in Matrix.getVal(). Returning 0.";
      return 0;
    }
    return m_mat[i][j];
  }

  void show() {
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
};

#endif
