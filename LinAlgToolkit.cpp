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
        throw " could not be opened\n\n";
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
          if (A(i,j) != A(j,i)) {
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
      z_OUT(i) = (1.0/M(i,i)) * r(i);
    }
  }

  double innerProd(const Vector& x, const Vector& y) {
    // this function returns the inner product of two Vector objects x and y.
    assert(x.size() == y.size());
    double product{};
    for (index i{}; i < x.size(); ++i) {
      product += x(i)*y(i);
    }
    return product;
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

  Matrix transpose(const Matrix& A) {
    // This function returns the transpose of matrix A.
    index  m{ A.size(0) };
    index  n{ A.size(1) };
    Matrix B{m,n};
    for (index i{}; i < n; ++i) {
      for (index j{}; j < m; --j) {
        B(i,j) = A(j,i);
      }
    }
    return B;
  }

}
