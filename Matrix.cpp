// This file contains the implementations of the member functions for the
// Matrix class
#include "Constants.h"
#include "LinAlgToolkit.h"
#include "Matrix.h"
#include "Types.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void Matrix::resize(const index m, const index n) {
  // to reset the size of the matrix
  m_row = m;
  m_col = n;
  m_mat.resize(m_row,vec_t(m_col));
}

index Matrix::size(const int dim) const {
  // to get the size of the matrix
  assert(dim==0 || dim==1);
  if (dim == 0) return m_mat.size();
  else          return m_mat[0].size();
}

void Matrix::setCol(const index j, const Vector x) {
  // to set the jth column of equal to a vector
  assert(m_row == x.size());
  for (index i{}; i < m_row; ++i) {
    m_mat[i][j] = x(i);
  }
}

void Matrix::show() const {
  // to print the matrix to the screen
  using std::cout;
  cout << std::setprecision(consts::dispPrec);
  cout << '\n';
  for (index i{}; i < m_row; ++i) {
    cout << '\n';
    for (auto a_ij : m_mat[i]) {
      cout << std::setw(consts::dispPrec+2)
           << std::right << std::showpoint << a_ij << '\t';
    }
  }
  cout << '\n';
}

void Matrix::load(const std::string filename) {
  // this function reads the contents of "filename" and stores the data in
  // the matrix A_OUT.
  // we'll read from a file called "filename"
  std::ifstream inf{filename};
  // in case we cannot open the input file stream
  tlk::catchFileError(filename, !inf);
  // get dimensions of matrix from file
  index m{};
  index n{};
  std::string strVal;
  inf >> strVal;
  m = static_cast<index>(std::stoi(strVal));
  inf >> strVal;
  n = static_cast<index>(std::stoi(strVal));
  // check that size of A_OUT matches size of matrix in file
  if (m != m_mat.size() || n != m_mat[0].size()) {
    std::cout << "\nError: Size of matrix A_OUT in readMatrix() does not "
              << "match size of matrix in " << filename << ".\n";
    std::exit(0);
  }
  // loop through file inserting values into matrix A
  for (index k{}; k < m*n; ++k) {
    // get value and store within A
    index i{ static_cast<index>(std::floor(k/n)) };
    index j{ k % n };
    inf >> strVal;
    m_mat[i][j] = std::stod(strVal);
  }
}

void Matrix::write(const std::string filename) const {
  // this function writes the contents of a matrix A into "filename"
  //create a file
  std::ofstream outf{filename};
  // in case we cannot open the output file stream
  if (!outf) {
    std::cerr << "\nError: soln.txt could not be opened for writing\n\n";
    return;
  }
  // get size of matrix
  index m{    m_mat.size() };
  index n{ m_mat[0].size() };
  // write into the file
  outf << m << ' ';
  outf << n << '\n';
  std::setprecision(consts::writePrec);
  for (index i{}; i < m; ++i) {
    for (index j{}; j < n; ++j) {
      outf << std::setw(consts::writePrec+2)
           << std::left << std::showpoint << m_mat[i][j] << ' ';
    }
    outf << '\n';
  }
  // when outf goes out of scope, the ofstream destructor will close the file
}

//*****************************************************************************/
// here are member functions for operator overloading
//*****************************************************************************/

double& Matrix::operator()(const index i, const index j) {
  // here we allow subscripting on our Matrix object
  assert(i < m_row && j < m_col);
  return m_mat[i][j];
}

double  Matrix::operator()(const index i, const index j) const {
  // here we allow subscripting on *constant* Matrix objects
  assert(i < m_row && j < m_col);
  return m_mat[i][j];
}

// Matrix& Matrix::operator=(const Matrix& A) {
//   // here we create an assignment operator
//   if (this == &A) return *this; // self-assignment guard
//
//   m_row = A.m_row;
//   m_col = A.m_col;
//   m_mat = A.m_mat;
//
//   return *this;
// }

//*****************************************************************************/
// here are regular functions for operator overloading
//*****************************************************************************/

Matrix operator-(const Matrix& A) {
  // here we negate all the elements of a Matrix A
  Matrix minus_A{A.size(0),A.size(1)};
  for (index i{}; i < A.size(0); ++i) {
    for (index j{}; j < A.size(1); ++j) {
        minus_A(i,j) = -A(i,j);
    }
  }
  return minus_A;
}

Matrix operator+(const Matrix& A, const Matrix& B) {
  // here we take the componentwise sum of two Matrices A and B
  assert(A.size(0) == B.size(0) && A.size(1) == B.size(1));
  Matrix C{A.size(0),A.size(1)};
  for (index i{}; i < A.size(0); ++i) {
    for (index j{}; j < A.size(1); ++j) {
        C(i,j) = A(i,j) + B(i,j);
    }
  }
  return C;
}

Matrix operator-(const Matrix& A, const Matrix& B) {
  // here we take the componentwise difference of two Matrices A and B
  assert(A.size(0) == B.size(0) && A.size(1) == B.size(1));
  return A + (-B); // we can now use the Vector operators defined above
}

Matrix operator*(const double scalar, const Matrix& A) {
  // here we scale a matrix
  Matrix scaled_A{A.size(0), A.size(1)};
  for (index i{}; i < A.size(0); ++i) {
    for (index j{}; j < A.size(1); ++j) {
      scaled_A(i,j) = scalar*A(i,j);
    }
  }
  return scaled_A;
}

Matrix operator*(const Matrix& A, const double scalar) {
  // we want scaling to be symmetric
  return scalar*A;
}

Vector operator*(const Matrix& A, const Vector& x) {
  // here we perform matrix-vector multiplication
  assert(A.size(1) == x.size());
  Vector y{x.size()};
  for (index i{}; i < x.size(); ++i) {
    for (index j{}; j < A.size(1); ++j) {
      y(i) = A(i,j)*x(j);
    }
  }
  return y;
}

Matrix operator*(const Matrix& A, const Matrix& B) {
  // here we perform matrix-matrix multiplication
  assert(A.size(1) == B.size(0));
  Matrix C{A.size(0), B.size(1)};
  for (index i{}; i < A.size(0); ++i) {
    for (index j{}; j < B.size(1); ++j) {
      for (index k{}; k < B.size(1); ++k) {
        C(i,j) += A(i,k)*B(k,j);
      }
    }
  }
  return C;
}
