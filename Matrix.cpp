// This file contains the implementations of the member functions for the
// Matrix class as well as its operator overloading.
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

//*****************************************************************************/
// Here are Matrix member functions
//*****************************************************************************/

void Matrix::resize(const index m, const index n) {
  // to reset the size of the matrix
  m_row = m;
  m_col = n;
  // m_mat.resize(m,vec_t(n)); // apparently this method doesn't work for resizing
  m_mat.resize(m);
  for (index i{}; i < m; ++i) {
    m_mat[i].resize(n);
  }
}

index Matrix::size(const int dim) const {
  // to get the size of the matrix
  assert(dim==0 || dim==1);
  if (dim == 0) return m_mat.size();
  else          return m_mat[0].size();
}

void Matrix::setCol(const index j, const Vector x) {
  // to set the jth column equal to a vector
  assert(m_row == x.size());
  for (index i{}; i < m_row; ++i) {
    m_mat[i][j] = x(i);
  }
}

void Matrix::setRow(const index i, const Vector x) {
  // to set the ith row equal to a vector
  assert(m_col == x.size());
  for (index j{}; j < m_col; ++j) {
    m_mat[i][j] = x(j);
  }
}

Vector Matrix::getCol(const index j) {
  // to extract the jth column as a Vector
  Vector x{m_row};
  for (index i{}; i < m_row; ++i) {
    x(i) = m_mat[i][j];
  }
  return x;
}

Vector Matrix::getRow(const index i) {
  // to extract the ith row as a Vector
  Vector x{m_col};
  for (index j{}; j < m_col; ++j) {
    x(j) = m_mat[i][j];
  }
  return x;
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
  // to reads the contents of "filename" and store the data in m_mat[][]
  std::ifstream inf{filename};
  tlk::catchFileError(filename, !inf);
  // get size of matrix from file
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
  // to write the contents of m_mat[][] into "filename"
  std::ofstream outf{filename};
  tlk::catchFileError(filename, !outf);
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
}

//*****************************************************************************/
// Here are Matrix member functions for operator overloading
//*****************************************************************************/

double& Matrix::operator()(const index i, const index j) {
  // to allow subscripting on our Matrix object
  assert(i < m_row && j < m_col);
  return m_mat[i][j];
}

double  Matrix::operator()(const index i, const index j) const {
  // to allow subscripting on *constant* Matrix objects
  assert(i < m_row && j < m_col);
  return m_mat[i][j];
}

// Matrix& Matrix::operator=(const Matrix& A) {
//   // to create an assignment operator
//   if (this == &A) return *this; // self-assignment guard
//
//   m_row = A.m_row;
//   m_col = A.m_col;
//   m_mat = A.m_mat;
//
//   return *this;
// }

//*****************************************************************************/
// Here are normal functions for operator overloading
//*****************************************************************************/

Matrix operator-(const Matrix& A) {
  // to negate all the elements of a Matrix A
  Matrix minus_A{A.size(0),A.size(1)};
  for (index i{}; i < A.size(0); ++i) {
    for (index j{}; j < A.size(1); ++j) {
        minus_A(i,j) = -A(i,j);
    }
  }
  return minus_A;
}

Matrix operator+(const Matrix& A, const Matrix& B) {
  // to take the componentwise sum of two Matrices A and B
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
  // to take the componentwise difference of two Matrices A and B
  assert(A.size(0) == B.size(0) && A.size(1) == B.size(1));
  return A + (-B); // we can now use the Vector operators defined above
}

Matrix operator*(const double scalar, const Matrix& A) {
  // to scale a matrix
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
  // to perform matrix-vector multiplication
  assert(A.size(1) == x.size());
  Vector y{A.size(0)};
  for (index i{}; i < A.size(0); ++i) {
    for (index j{}; j < A.size(1); ++j) {
      y(i) += A(i,j)*x(j);
    }
  }
  return y;
}

Matrix operator*(const Matrix& A, const Matrix& B) {
  // to perform matrix-matrix multiplication
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

std::ostream& operator<<(std::ostream& out, const Matrix& A) {
  // to output Matrix objects easily
  index m{A.size(0)};
  index n{A.size(1)};
  out << "\n" << std::setprecision(consts::dispPrec);

  // show all or partial Matrix depending on its size
  index max_size{ std::max(m, n) };
  index min_size{ std::min(m, n) };

  // if small enough, show as usual
  if (max_size <= consts::max_out_size) {
    for (index i{}; i < m; ++i) {
      out << '\n';
      for (index j{}; j < n; ++j) {
        out << std::setw(consts::dispPrec+2) << std::right
                                      << std::showpoint << A(i,j) << '\t';
      }
    }
  }

  // if not, we loop only over max_out_size-many elements and print by case
  else if (min_size > consts::max_out_size) {
    for (index i{}; i < consts::max_out_size; ++i) {
      out << '\n';
      for (index j{}; j < consts::max_out_size; ++j) {

        // Case 1: upper left corner => print as usual
        if (i < consts::max_out_size/2 && j < consts::max_out_size/2) {
          out << std::setw(consts::dispPrec+2) << std::right
              << std::showpoint << A(i,j) << '\t';
        }
        // Case 2: upper right corner => print as usual, but include "..."
        //         between left and right blocks
        else if (i < consts::max_out_size/2 && j >= consts::max_out_size/2) {
          if (j == consts::max_out_size/2) out << ". . .\t";
          out << std::setw(consts::dispPrec+2) << std::right << std::showpoint
              << A( i, A.size(1) -consts::max_out_size +j ) << '\t';
        }
        // Case 3: lower left corner => print as usual but include vertical
        //         "..." between upper and lower blocks
        else if (i >= consts::max_out_size/2 && j < consts::max_out_size/2) {
          if (i == consts::max_out_size/2 && j == 0) {
            for (index l{}; l < 3; ++l) {
              out << '\t';
              for (index k{}; k < consts::max_out_size; ++k) {
                if (k == consts::max_out_size/2 - 1) {
                  out << ".\t";
                }
                else if (k == consts::max_out_size/2) {
                  if      (l == 0) out << ".    \t\t";
                  else if (l == 1) out << "  .  \t\t";
                  else             out << "    .\t\t";
                }
                else out << ".\t\t";
              }
              out << '\n';
            }
          }
          out << std::setw(consts::dispPrec+2) << std::right << std::showpoint
              << A( A.size(0) -consts::max_out_size +i, j ) << '\t';
        }
        // Case 4: lower right corner => print as usual
        else {
          if (j == consts::max_out_size/2) out << ". . .\t";
          out << std::setw(consts::dispPrec+2) << std::right << std::showpoint
              << A( A.size(0) -consts::max_out_size +i,
                    A.size(1) -consts::max_out_size +j ) << '\t';
        }
      }
    }
  }
  return out << "\n"; // return std::ostream
}
