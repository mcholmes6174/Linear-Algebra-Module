/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains the entire contents of the namespace "tlk::". The
 * functions within this namespace are written to be used when doing computation
 * with matrices and vectors, particularly in the context of numerical linear
 * algebra.
 *
 * Functions that would otherwise return Vector or Matrix objects utilize
 * non-const reference parameters for performance reasons. These parameters are
 * labeled with the suffix _OUT. Generally, all other functions return by value
 * unless void.
 *
 */
#include "Constants.h" // consts:: user-defined namespace
#include "Matrix.h"
#include "Types.h"
#include "Vector.h"
#include <cassert>     // for assert() function calls
#include <cmath>       // std::pow(), std::sqrt(), std::abs(), etc.
#include <fstream>     // to do file I/O
#include <iomanip>     // to override the default precision shown by std::cout
#include <iostream>
#include <string>

namespace tlk {

  void catchFileError(const std::string filename, const bool open_failed) {
    // exception handling function for input file stream.
    try {
      if (open_failed) {
        throw " could not be opened for reading\n\n";
      }
      else {
        return;
      }
    }
    catch (const char* exception) {
      std::cerr << "\nError: " << filename << exception;
      std::exit(0);
    }
  }

  void catchSingular(const double matrix_entry) {
    // exception handling funciton for testing for matrix singularity.
    try {
      if (std::abs(matrix_entry) < consts::epsilon) {
        throw "Matrix is singular! Exiting the program.";
      }
      else {
        return;
      }
    }
    catch (const char* exception) {
      std::cerr << "\nError: " << exception << '\n' << '\n';
      std::exit(0);
    }
  }

  void checkSymm(const Matrix& A) {
    // this function checks whether a given square matrix is symmetric, and
    // throws an error if it is not symmetric.
    for (index i{}; i < A.size(0); ++i) {
      for (index j{}; j < i; ++j) {
        try {
          if (A.get(i,j) != A.get(j,i)) {
            throw "Matrix is not symmetric! C-G method will not be successful.";
          }
        }
        catch (const char* exception) {
          std::cerr << "\nError: " << exception << '\n' << '\n';
          std::exit(0);
        }
      }
    }
  }

  void diagPreCond(const Matrix& M, const Vector& r, Vector& z_OUT) {
    // this function performs diagonal preconditioning for C-G method
    // via the multiplication z = M^(-1)*r, where M=diag(A) is the diagonal
    // preconditoner matrix. The original matrix A should be passed as the first
    // argument.
    assert(r.size() == z_OUT.size() && r.size() == M.size(0));
    for (index i{}; i < r.size(); ++i) {
      z_OUT.set(i) = (1.0/M.get(i,i)) * r.get(i);
    }
  }

  void getError(const Matrix& A, const Vector& x, const Vector& b,
                                                        Vector& err_OUT) {
    // this function computes the error in the solution x to the system Ax=b
    // once a solution has been found. It simply returns err = Ax-b.
    // forward declaration necessary if we are to keep alphabetical order
    void matVecMul(const Matrix&, const Vector&, Vector& y_OUT);
    matVecMul(A,x,err_OUT);       // multiply A*x
    for (index i{}; i < x.size(); ++i) {
      err_OUT.set(i) -= b.get(i); // subtract b
    }
  }

  double innerProd(const Vector& x, const Vector& y) {
    // this function returns the inner product of two Vector objects x and y.
    assert(x.size() == y.size());
    double product{};
    for (index i{}; i < x.size(); ++i) {
      product += x.get(i)*y.get(i);
    }
    return product;
  }

  void matMul(const Matrix& A, const Matrix& B, Matrix& C_OUT) {
    // this function performs the matrix multiplication C=AB.
    assert(A.size(1) == B.size(0));
    assert(A.size(0) == C_OUT.size(0));
    assert(B.size(1) == C_OUT.size(1));
    for (index i{}; i < A.size(0); ++i) {
      for (index j{}; j < B.size(1); ++j) {
        for (index k{}; k < B.size(0); ++k) {
          C_OUT.set(i,j) += A.get(i,k) * B.get(k,j);
        }
      }
    }
  }

  Vector makeVecCopy(const Vector& x, const double scalar=1.0) {
    // this function returns the vector x after optionally being scaled.
    index m{ x.size() };
    Vector x_copy{m};
    for (index i{}; i < m; ++i) {
      x_copy.set(i) = scalar*x.get(i);
    }
    // we return by value here so that a copy is made automatically
    return x_copy;
  }

  void matVecMul(const Matrix& A, const Vector& x, Vector& y_OUT) {
    // this function performs the matrix-vector multiplication y=Ax.
    assert(A.size(0) == y_OUT.size() && A.size(1) == x.size());
    for (index i{}; i < A.size(0); ++i) {
      for (index j{}; j < A.size(1); ++j) {
        y_OUT.set(i) += A.get(i,j) * x.get(j);
      }
    }
  }

  void nthColumn(const std::string action, Matrix& A,
                                           Vector& x, const index n) {
    // This function either returns the nth column of a matrix as a vector or
    // stores a vector in the nth column of the given matrix based on the value
    // of the string "action".
    assert(A.size(0) == x.size());
    if (action == "grab") {
      for (index i{}; i < x.size(); ++i) {
        x.set(i) = A.get(i,n);
      }
    }
    else if (action == "give") {
      for (index i{}; i < x.size(); ++i) {
        A.set(i,n) = x.get(i);
      }
    }
  }

  index readNthVal(const std::string filename, const int n) {
    // this function reads the nth value from the file "filename" and returns
    // it as type std::index. This function was created in order to read the
    // dimensions of a matrix in from file.
    std::ifstream inf{filename};
    catchFileError(filename,!inf);
    std::string strVal;
    for (int i{}; i <= n; ++i) {
      inf >> strVal;
    }
    return static_cast<index>( std::stoi(strVal) );
  }

  void swapRows(Matrix& A_OUT, Vector& b_OUT, const index r1, const index r2) {
    // This function swaps row1 and row2 of both the matrix A
    // and vector b in the system Ax=b.
    assert(A_OUT.size(0) == b_OUT.size());
    double placeholder{};
    // swap rows in matrix
    for (index k{}; k < A_OUT.size(1); ++k) {
      placeholder     = A_OUT.get(r1,k);
      A_OUT.set(r1,k) = A_OUT.get(r2,k);
      A_OUT.set(r2,k) = placeholder;
    }
    // swap rows in vector
    placeholder   = b_OUT.get(r1);
    b_OUT.set(r1) = b_OUT.get(r2);
    b_OUT.set(r2) = placeholder;
  }

  Matrix transpose(const Matrix& A) {
    // This function returns the transpose of matrix A.
    index  m{ A.size(0) };
    index  n{ A.size(1) };
    Matrix B{m,n};
    for (index i{}; i < n; ++i) {
      for (index j{}; j < m; --j) {
        B.set(i,j) = A.get(j,i);
      }
    }
    return B;
  }

  Vector updateVector(const Vector& vec1, const double scalar,
                                          const Vector& vec2) {
    // This function returns the linear update x = vec1 + scalar*vec2.
    // Here we choose to return by value (rather than pass x_OUT by reference)
    // in order to increase readability.
    assert(vec1.size() == vec2.size());
    index  m{vec1.size()};
    Vector x{m};
    for (index i{}; i < m; ++i) {
      x.set(i) = vec1.get(i) + scalar*vec2.get(i);
    }
    return x;
  }

}
