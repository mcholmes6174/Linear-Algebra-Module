/* Written by Matthew Holmes
 * November 2021
 *
 * This program contains the entire contents of the namespace "tools::". The
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

namespace tools {

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
    for (int i{}; i < std::size(A); ++i) {
      for (int j{}; j < i; ++j) {
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
    int m{std::size(M)};
    vec_t z(m);
    for (int i{}; i < m; ++i) {
      z[i] = (1.0/M[i][i])*r[i];
    }
    return z;
  }

  vec_t getError(const mat_t& A, const vec_t& x, const vec_t& b) {
    // this function computes the error in the solution x to the system Ax=b
    // once a solution has been found. It simply returns err = Ax-b.
    // need to forward declare matVecMul here
    vec_t matVecMul(const mat_t&, const vec_t&);
    int m{std::size(A)};
    vec_t err(m);
    err = matVecMul(A,x);  // multiply A*x
    for (int i{}; i < m; ++i) {
      err[i] -= b[i];      // subtract b
    }
  }

  double getNorm(const vec_t& x) {
    // this function returns the Euclidean norm of a vector x.
    double norm{};
    // for (int i{}; i < std::size(x); ++i) { // need to size size() bug!!!!
    for (int i{}; i < 4; ++i) {
      norm += std::pow(x[i],2.0);
    }
    return std::sqrt(norm);
  }

  double innerProd(const vec_t& x, const vec_t& y) {
    // this function returns the inner product of two vectors x and y.
    if (std::size(x) != std::size(y)) {
      std::cerr << "Error: Vectors must have equal length to compute <x,y>";
      std::exit(0);
    }
    double result{};
    for (int i{}; i < std::size(x); ++i) {
      result += x[i]*y[i];
    }
    return result;
  }

  mat_t matMul(const mat_t& A, const mat_t& B) {
    // this function performs the matrix multiplication C=AB.
    int m{std::size(A)};
    int n{std::size(A[0])};
    mat_t C(m,vec_t(n));
    for (int i{}; i < m; ++i) {
      for (int j{}; j < n; ++j) {
        for (int k{}; k < n; ++k) {
          C[i][j] += A[i][k] * B[k][j];
        }
      }
    }
    return C;
  }

  vec_t makeVecCopy(const vec_t& x) {
    // this function returns a copy of the vector x.
    int m{std::size(x)};
    vec_t x_copy(m);
    for (int i{}; i < m; ++i) {
      x_copy[i] = x[i];
    }
    return x_copy;
  }

  vec_t matVecMul(const mat_t& A, const vec_t& x) {
    // this function performs the matrix-vector multiplication y=Ax.
    int m{std::size(A)};
    vec_t y(m);
    for (int i{}; i < m; ++i) {
      for (int j{}; j < std::size(A[0]); ++j) {
        y[i] += A[i][j] * x[j];
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
    int m{};
    int n{};
    std::string strVal;
    inf >> strVal;
    m = std::stoi(strVal);
    inf >> strVal;
    n = std::stoi(strVal);
    // declare m by n matrix
    mat_t A(m,vec_t(n));
    // loop through file inserting values into matrix A
    for (int k{}; k < m*n; ++k) {
      // get value and store within A
      int i{ static_cast<int>(std::floor(k/m)) };
      int j{ k % n };
      inf >> strVal;
      A[i][j] = std::stod(strVal);
    }
    return A;
  }

  int readNthVal(const std::string filename, const int n) {
    // this function reads the nth value from the file "filename" and
    // returns it as type integer
    std::ifstream inf{filename};
    catchFileError(filename,!inf);
    std::string strVal;
    for (int i{}; i <= n; ++i) {
      inf >> strVal;
    }
    return std::stoi(strVal);
  }

  vec_t readVector(const std::string filename) {
    // this function reads the contents of "filename" and stores the data in
    // the vector x, which is returned on exit.
    // we'll read from a file called "filename"
    std::ifstream inf{filename};
    // in case we cannot open the input file stream
    catchFileError(filename, !inf);
    // get length of vector from file
    int m{};
    std::string strVal;
    inf >> strVal;
    m = std::stoi(strVal);
    // declare vector
    vec_t x(m);
    // loop through file inserting values into vector
    for (int i{}; i < m; ++i) {
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
    for (int i{}; i < std::size(A); ++i) {
      cout << '\n';
      for (int j{}; j < std::size(A[0]); ++j) {
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
    for (int i{}; i < std::size(x); ++i) {
      cout << '\n' << std::setw(dispPrec+2) << std::right
                                            << std::showpoint << x[i];
    }
    cout << '\n';
  }

  void swapRows(mat_t& A, vec_t& b, const int r1, const int r2) {
    // This function swaps row1 and row2 of both the matrix A and vector b in
    // the system Ax=b
    double placeholder{};
    // swap rows in matrix
    for (int k{}; k < std::size(A[0]); ++k) {
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
    int m{std::size(A)};
    int n{std::size(A[0])};
    // write into the file
    outf << m << '\n';
    outf << n << '\n';
    std::setprecision(writePrec);
    for (int i{}; i < m; ++i) {
      for (int j{}; j < n; ++j) {
        outf << std::setw(writePrec+2) << std::right
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
    // write into the file
    std::setprecision(writePrec);
    for (int i{}; i < std::size(x); ++i) {
      outf << std::setw(writePrec+2) << std::right
                                     << std::showpoint << x[i] << '\n';
    }
    // when outf goes out of scope, the ofstream destructor will close the file
  }

}
