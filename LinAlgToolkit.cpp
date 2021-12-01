/* Written by Matthew Holmes
 * November 2021
 *
 * This program contains the entire contents of the namespace "tlk::". The
 * functions within this namespace are written to be used when doing computation
 * with matrices and vectors, particularly in the context of numerical lienar
 * algebra.
 *
 */
#include "Constants.h" // for consts::XXX user-defined namespace
#include "Types.h"     // for type alias mat_t and vec_t
#include <cmath>       // for std::pow() and std::sqrt() and std::abs()
#include <fstream>     // to do file I/O
#include <iomanip>     // to override the default precision shown by std::cout
#include <iostream>
#include <string>
#include <vector>

namespace tlk {

  using std::size_t; // unsigned int type for std::vector sizes and indexing

  int dispPrec{8};   // precision with which to display double values on screen
  int writePrec{16}; // precision with which to write double values to .dat file

  void catchFileError(const std::string filename, const bool open_failed) {
    // Exception handling function for input file stream
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
    // Exception handling funciton for testing for matrix singularity.
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

  void checkSymm(const mat_t& A) {
    // this function checks whether a given square matrix is symmetric, and
    // throws an error if it is not symmetric.
    for (size_t i{}; i < A.size(); ++i) {
      for (size_t j{}; j < i; ++j) {
        try {
          if (A[i][j] != A[j][i]) {
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

  vec_t diagPreCond(const mat_t& M, const vec_t& r) {
    // this function performs diagonal preconditioning for C-G method
    // via the multiplication z = M^(-1)*r, where M=diag(A) is the diagonal
    // preconditoner matrix. The original matrix A should be passed as the first
    // argument.
    size_t m{M.size()};
    vec_t  z(m);
    for (size_t i{}; auto r_i : r) {
      z[i] = (1.0/M[i][i]) * r_i;
      ++i;
    }
    return z;
  }

  vec_t getError(const mat_t& A, const vec_t& x, const vec_t& b) {
    // this function computes the error in the solution x to the system Ax=b
    // once a solution has been found. It simply returns err = Ax-b.
    size_t m{A.size()};
    vec_t  err(m);
    // need to forward declare matVecMul
    vec_t matVecMul(const mat_t&, const vec_t&);
    err = matVecMul(A,x);  // multiply A*x
    for (size_t i{}; auto b_i : b) {
      err[i] -= b_i;       // subtract b
      ++i;
    }
    return err;
  }

  double getNorm(const vec_t& x) {
    // this function returns the Euclidean norm of a vector x.
    double norm{};
    // here we use a for-each loop
    for (auto x_i : x) {
      norm += std::pow(x_i,2.0);
    }
    return std::sqrt(norm);
  }

  double innerProd(const vec_t& x, const vec_t& y) {
    // this function returns the inner product of two vectors x and y.
    if (x.size() != y.size()) {
      std::cerr << "Error: Vectors must have equal length to compute <x,y>";
      std::exit(0);
    }
    double result{};
    for (size_t i{}; auto x_i : x) {
      result += x_i*y[i];
      ++i;
    }
    return result;
  }

  mat_t matMul(const mat_t& A, const mat_t& B) {
    // this function performs the matrix multiplication C=AB.
    size_t m{A.size()};
    size_t n{A[0].size()};
    mat_t C(m,vec_t(n));
    for (size_t i{}; i < m; ++i) {
      for (size_t j{}; j < n; ++j) {
        for (size_t k{}; k < n; ++k) {
          C[i][j] += A[i][k] * B[k][j];
        }
      }
    }
    return C;
  }

  vec_t makeVecCopy(const vec_t& x) {
    // this function returns a copy of the vector x.
    size_t m{x.size()};
    vec_t  x_copy(m);
    for (size_t i{}; auto x_i : x) {
      x_copy[i] = x_i;
      ++i;
    }
    return x_copy;
  }

  vec_t matVecMul(const mat_t& A, const vec_t& x) {
    // this function performs the matrix-vector multiplication y=Ax.
    size_t m{A.size()};
    vec_t  y(m);
    for (size_t i{}; i < m; ++i) {
      for (size_t j{}; auto x_j : x) {
        y[i] += A[i][j] * x_j;
        ++j;
      }
    }
    return y;
  }

  mat_t readMatrix(const std::string filename) {
    // this function reads the contents of "filename" and stores the data in
    // the matrix A, which is returned on exit.
    // we'll read from a file called "filename"
    std::ifstream inf{filename};
    // in case we cannot open the input file stream
    catchFileError(filename, !inf);
    // get dimensions of matrix from file
    size_t m{};
    size_t n{};
    std::string strVal;
    inf >> strVal;
    m = static_cast<size_t>(std::stoi(strVal));
    inf >> strVal;
    n = static_cast<size_t>(std::stoi(strVal));
    // declare m by n matrix
    mat_t A(m,vec_t(n));
    // loop through file inserting values into matrix A
    for (size_t k{}; k < m*n; ++k) {
      // get value and store within A
      size_t i{ static_cast<size_t>(std::floor(k/m)) };
      size_t j{ k % n };
      inf >> strVal;
      A[i][j] = std::stod(strVal);
    }
    return A;
  }

  size_t readNthVal(const std::string filename, const int n) {
    // this function reads the nth value from the file "filename" and
    // returns it as type integer
    std::ifstream inf{filename};
    catchFileError(filename,!inf);
    std::string strVal;
    for (int i{}; i <= n; ++i) {
      inf >> strVal;
    }
    return static_cast<size_t>(std::stoi(strVal));
  }

  vec_t readVector(const std::string filename) {
    // this function reads the contents of "filename" and stores the data in
    // the vector x, which is returned on exit.
    // we'll read from a file called "filename"
    std::ifstream inf{filename};
    // in case we cannot open the input file stream
    catchFileError(filename, !inf);
    // get length of vector from file
    size_t m{};
    std::string strVal;
    inf >> strVal;
    m = static_cast<size_t>(std::stoi(strVal));
    // declare vector
    vec_t x(m);
    // loop through file inserting values into vector
    for (size_t i{}; i < m; ++i) {
      // get value and store within A
      inf >> strVal;
      x[i] = std::stod(strVal);
    }
    return x;
  }

  void showMatrix(const mat_t& A) {
    // this function prints a matrix to the screen.
    using std::cout;
    cout << std::setprecision(dispPrec);
    for (size_t i{}; i < A.size(); ++i) {
      cout << '\n';
      for (size_t j{}; j < A[0].size(); ++j) {
        cout << std::setw(dispPrec+2) << std::right
                                      << std::showpoint << A[i][j] << '\t';
      }
    }
    cout << '\n';
  }

  void showVector(const vec_t& x) {
    // this function prints a vector to the screen.
    using std::cout;
    cout << std::setprecision(dispPrec);
    for (auto x_i : x) {
      cout << '\n' << std::setw(dispPrec+2) << std::right
                                            << std::showpoint << x_i;
    }
    cout << '\n';
  }

  void swapRows(mat_t& A, vec_t& b, const size_t r1, const size_t r2) {
    // This function swaps row1 and row2 of both the matrix A and vector b in
    // the system Ax=b
    double placeholder{};
    // swap rows in matrix
    for (size_t k{}; k < A[0].size(); ++k) {
      placeholder = A[r1][k];
      A[r1][k]    = A[r2][k];
      A[r2][k]    = placeholder;
    }
    // swap rows in vector
    placeholder = b[r1];
    b[r1]       = b[r2];
    b[r2]       = placeholder;
  }

  void writeMatrix(const std::string filename, const mat_t& A) {
    // this function writes the contents of a matrix A into a file "filename"
    // create a file
    std::ofstream outf{filename};
    // in case we cannot open the output file stream
    if (!outf) {
      std::cerr << "\nError: soln.txt could not be opened for writing\n\n";
      return;
    }
    // get size of matrix
    size_t m{A.size()};
    size_t n{A[0].size()};
    // write into the file
    outf << m << ' ';
    outf << n << '\n';
    std::setprecision(writePrec);
    for (size_t i{}; i < m; ++i) {
      for (size_t j{}; j < n; ++j) {
        outf << std::setw(writePrec+2) << std::left
                                       << std::showpoint << A[i][j] << ' ';
      }
      outf << '\n';
    }
    // when outf goes out of scope, the ofstream destructor will close the file
  }

  void writeVector(const std::string filename, const vec_t& x) {
    // this function writes the contents of a vector x into a file "filename"
    // create a file
    std::ofstream outf{filename};
    // in case we cannot open the output file stream
    if (!outf) {
      std::cerr << "\nError: soln.txt could not be opened for writing\n\n";
      return;
    }
    // get length of vector and write into file
    size_t m{x.size()};
    outf << m << '\n';
    // write into the file
    std::setprecision(writePrec);
    for (auto x_i : x) {
      outf << std::setw(writePrec+2) << std::left
                                     << std::showpoint << x_i << '\n';
    }
    // when outf goes out of scope, the ofstream destructor will close the file
  }

}
