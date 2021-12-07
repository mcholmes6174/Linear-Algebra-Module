/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains the entire contents of the namespace "tlk::". The
 * functions within this namespace are written to be used when doing computation
 * with matrices and vectors, particularly in the context of numerical linear
 * algebra.
 *
 * Functions that would otherwise return vec_t or mat_t types utilize non-const
 * reference parameters for performance reasons. These parameters are labeled
 * with the suffix _OUT. All other functions return by value unless void.
 *
 */
#include "Constants.h" // consts:: user-defined namespace
#include "Types.h"
#include <cmath>       // std::pow() and std::sqrt() and std::abs()
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

  void diagPreCond(const mat_t& M, const vec_t& r, vec_t& z_OUT) {
    // this function performs diagonal preconditioning for C-G method
    // via the multiplication z = M^(-1)*r, where M=diag(A) is the diagonal
    // preconditoner matrix. The original matrix A should be passed as the first
    // argument.
    for (size_t i{}; auto r_i : r) {
      z_OUT[i] = (1.0/M[i][i]) * r_i;
      ++i;
    }
  }

  void getError(const mat_t& A, const vec_t& x, const vec_t& b, vec_t& err_OUT) {
    // this function computes the error in the solution x to the system Ax=b
    // once a solution has been found. It simply returns err = Ax-b.
    void matVecMul(const mat_t&, const vec_t&, vec_t& y_OUT); // forward declare
    matVecMul(A,x,err_OUT);  // multiply A*x
    for (size_t i{}; auto b_i : b) {
      err_OUT[i] -= b_i;     // subtract b
      ++i;
    }
  }

  double getNorm(const vec_t& x) {
    // this function returns the Euclidean norm of a vector x.
    double norm{};
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
    double product{};
    for (size_t i{}; auto x_i : x) {
      product += x_i*y[i];
      ++i;
    }
    return product;
  }

  void matMul(const mat_t& A, const mat_t& B, mat_t C_OUT) {
    // this function performs the matrix multiplication C=AB.
    if (A[0].size() != B.size()) {
      std::cerr << "Error: Dimensions must be consistent for matMul()";
      std::exit(0);
    }
    for (size_t i{}; i < A.size(); ++i) {
      for (size_t j{}; j < B.size(); ++j) {
        for (size_t k{}; k < B.size(); ++k) {
          C_OUT[i][j] += A[i][k] * B[k][j];
        }
      }
    }
  }

  vec_t makeVecCopy(const vec_t& x, const double scalar=1.0) {
    // this function returns the vector x after optionally being scaled.
    vec_t x_OUT(x.size());
    for (size_t i{}; auto x_i : x) {
      x_OUT[i] = scalar*x_i;
      ++i;
    }
    // we return by value here so that a copy is made automatically
    return x_OUT;
  }

  void matVecMul(const mat_t& A, const vec_t& x, vec_t& y_OUT) {
    // this function performs the matrix-vector multiplication y=Ax.
    for (size_t i{}; i < A.size(); ++i) {
      for (size_t j{}; auto x_j : x) {
        y_OUT[i] += A[i][j] * x_j;
        ++j;
      }
    }
  }

  void normalizeVec(vec_t& x_OUT) {
    double norm{ getNorm(x_OUT) };
    for (size_t i{}; i < x_OUT.size(); ++i) {
      x_OUT[i] /= norm;
    }
  }

  void nthColumn(const std::string action, mat_t& A, vec_t& x, const size_t n) {
    // This function either returns the nth column of a matrix as a vector or
    // stores a vector in the nth column of the given matrix based on the value
    // of the string "action".
    if (action == "grab") {
      for (size_t i{}; auto a_i : A) {
        x[i] = a_i[n];
        ++i;
      }
    }
    else if (action == "give") {
      for (size_t i{}; auto x_i : x) {
        A[i][n] = x_i;
        ++i;
      }
    }
  }

  void readMatrix(const std::string filename, mat_t& A_OUT) {
    // this function reads the contents of "filename" and stores the data in
    // the matrix A_OUT.
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
    // check that size of A_OUT matches size of matrix in file
    if ( (m != A_OUT.size()) || (n != A_OUT[0].size()) ) {
      std::cout << "\nError: Size of matrix A_OUT in readMatrix() does not "
                << "match size of matrix in " << filename << ".\n";
      std::exit(0);
    }
    // loop through file inserting values into matrix A
    for (size_t k{}; k < m*n; ++k) {
      // get value and store within A
      size_t i{ static_cast<size_t>(std::floor(k/n)) };
      size_t j{ k % n };
      inf >> strVal;
      A_OUT[i][j] = std::stod(strVal);
    }
  }

  size_t readNthVal(const std::string filename, const int n) {
    // this function reads the nth value from the file "filename" and returns
    // it as type std::size_t. This function was created in order to read the
    // dimensions of a matrix in from file.
    std::ifstream inf{filename};
    catchFileError(filename,!inf);
    std::string strVal;
    for (int i{}; i <= n; ++i) {
      inf >> strVal;
    }
    return static_cast<size_t>( std::stoi(strVal) );
  }

  void readVector(const std::string filename, vec_t& x_OUT) {
    // this function reads the contents of "filename" and stores the data in
    // the vector x_OUT.
    // we'll read from a file called "filename"
    std::ifstream inf{filename};
    // in case we cannot open the input file stream
    catchFileError(filename, !inf);
    // get length of vector from file
    size_t m{};
    std::string strVal;
    inf >> strVal;
    m = static_cast<size_t>(std::stoi(strVal));
    // check that size of x_OUT matches size of vector in file
    if (m != x_OUT.size()) {
      std::cout << "\nError: Size of vector x_OUT in readVector() does not "
                << "match size of vector in " << filename << ".\n";
      std::exit(0);
    }
    // loop through file inserting values into vector
    for (size_t i{}; i < m; ++i) {
      // get value and store within A
      inf >> strVal;
      x_OUT[i] = std::stod(strVal);
    }
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

  void swapRows(mat_t& A_OUT, vec_t& b_OUT, const size_t r1, const size_t r2) {
    // This function swaps row1 and row2 of both the matrix A
    // and vector b in the system Ax=b.
    double placeholder{};
    // swap rows in matrix
    for (size_t k{}; k < A_OUT[0].size(); ++k) {
      placeholder  = A_OUT[r1][k];
      A_OUT[r1][k] = A_OUT[r2][k];
      A_OUT[r2][k] = placeholder;
    }
    // swap rows in vector
    placeholder = b_OUT[r1];
    b_OUT[r1]   = b_OUT[r2];
    b_OUT[r2]   = placeholder;
  }

  mat_t transpose(const mat_t& A) {
    // This function returns the transpose of matrix A.
    size_t m{    A.size() };
    size_t n{ A[0].size() };
    mat_t B(n,vec_t(m));
    for (size_t i{}; i < n; ++i) {
      for (size_t j{}; j < m; --j) {
        B[i][j] = A[j][i];
      }
    }
    return B;
  }

  vec_t updateVector(const vec_t& vec1, const double scalar, const vec_t& vec2) {
    // This function returns the linear update x = vec1 + scalar*vec2.
    // Here we choose to return by value (rather than pass x_OUT by reference)
    // in order to increase readability.
    size_t m{vec1.size()};
    vec_t  x(m);
    for (size_t i{}; auto vec2_i : vec2) {
      x[i] = vec1[i] + scalar*vec2_i;
      ++i;
    }
    return x;
  }

  void writeMatrix(const std::string filename, const mat_t& A) {
    // this function writes the contents of a matrix A into "filename"
    //create a file
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
    // this function writes the contents of a vector x into "filename"
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
