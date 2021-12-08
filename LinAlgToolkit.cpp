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
#include "Classes.h"
#include "Constants.h" // consts:: user-defined namespace
#include "Types.h"
#include <cmath>       // std::pow() and std::sqrt() and std::abs()
#include <fstream>     // to do file I/O
#include <iomanip>     // to override the default precision shown by std::cout
#include <iostream>
#include <string>
#include <vector>

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

  void checkSymm(const mat_t& A) {
    // this function checks whether a given square matrix is symmetric, and
    // throws an error if it is not symmetric.
    for (index i{}; i < A.size(); ++i) {
      for (index j{}; j < i; ++j) {
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
    for (index i{}; auto r_i : r) {
      z_OUT[i] = (1.0/M[i][i]) * r_i;
      ++i;
    }
  }

  void getError(const mat_t& A, const vec_t& x, const vec_t& b, vec_t& err_OUT) {
    // this function computes the error in the solution x to the system Ax=b
    // once a solution has been found. It simply returns err = Ax-b.
    void matVecMul(const mat_t&, const vec_t&, vec_t& y_OUT); // forward declare
    matVecMul(A,x,err_OUT);  // multiply A*x
    for (index i{}; auto b_i : b) {
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
    for (index i{}; auto x_i : x) {
      product += x_i*y[i];
      ++i;
    }
    return product;
  }

  double innerProd2(Vector& x, Vector& y) {
    // this function returns the inner product of two Vector objects x and y.
    if (x.getSize() != y.getSize()) {
      std::cerr << "Error: Vectors must have equal length to compute <x,y>";
      std::exit(0);
    }
    double product{};
    for (index i{}; i < x.getSize(); ++i) {
      product += x.getVal(i)*y.getVal(i);
    }
    return product;
  }

  void matMul(const mat_t& A, const mat_t& B, mat_t C_OUT) {
    // this function performs the matrix multiplication C=AB.
    if (A[0].size() != B.size()) {
      std::cerr << "Error: Dimensions must be consistent for matMul()";
      std::exit(0);
    }
    for (index i{}; i < A.size(); ++i) {
      for (index j{}; j < B.size(); ++j) {
        for (index k{}; k < B.size(); ++k) {
          C_OUT[i][j] += A[i][k] * B[k][j];
        }
      }
    }
  }

  vec_t makeVecCopy(const vec_t& x, const double scalar=1.0) {
    // this function returns the vector x after optionally being scaled.
    vec_t x_OUT(x.size());
    for (index i{}; auto x_i : x) {
      x_OUT[i] = scalar*x_i;
      ++i;
    }
    // we return by value here so that a copy is made automatically
    return x_OUT;
  }

  void matVecMul(const mat_t& A, const vec_t& x, vec_t& y_OUT) {
    // this function performs the matrix-vector multiplication y=Ax.
    for (index i{}; i < A.size(); ++i) {
      for (index j{}; auto x_j : x) {
        y_OUT[i] += A[i][j] * x_j;
        ++j;
      }
    }
  }

  void normalizeVec(vec_t& x_OUT) {
    double norm{ getNorm(x_OUT) };
    for (index i{}; i < x_OUT.size(); ++i) {
      x_OUT[i] /= norm;
    }
  }

  void nthColumn(const std::string action, mat_t& A, vec_t& x, const index n) {
    // This function either returns the nth column of a matrix as a vector or
    // stores a vector in the nth column of the given matrix based on the value
    // of the string "action".
    if (action == "grab") {
      for (index i{}; auto a_i : A) {
        x[i] = a_i[n];
        ++i;
      }
    }
    else if (action == "give") {
      for (index i{}; auto x_i : x) {
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
    index m{};
    index n{};
    std::string strVal;
    inf >> strVal;
    m = static_cast<index>(std::stoi(strVal));
    inf >> strVal;
    n = static_cast<index>(std::stoi(strVal));
    // check that size of A_OUT matches size of matrix in file
    if ( (m != A_OUT.size()) || (n != A_OUT[0].size()) ) {
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
      A_OUT[i][j] = std::stod(strVal);
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

  void readVector(const std::string filename, vec_t& x_OUT) {
    // this function reads the contents of "filename" and stores the data in
    // the vector x_OUT.
    // we'll read from a file called "filename"
    std::ifstream inf{filename};
    // in case we cannot open the input file stream
    catchFileError(filename, !inf);
    // get length of vector from file
    index m{};
    std::string strVal;
    inf >> strVal;
    m = static_cast<index>(std::stoi(strVal));
    // check that size of x_OUT matches size of vector in file
    if (m != x_OUT.size()) {
      std::cout << "\nError: Size of vector x_OUT in readVector() does not "
                << "match size of vector in " << filename << ".\n";
      std::exit(0);
    }
    // loop through file inserting values into vector
    for (index i{}; i < m; ++i) {
      // get value and store within A
      inf >> strVal;
      x_OUT[i] = std::stod(strVal);
    }
  }

  void showMatrix(const mat_t& A) {
    // this function prints a matrix to the screen.
    using std::cout;
    cout << std::setprecision(consts::dispPrec);
    for (index i{}; i < A.size(); ++i) {
      cout << '\n';
      for (auto a_ij : A[i]) {
        cout << std::setw(consts::dispPrec+2)
             << std::right << std::showpoint << a_ij << '\t';
      }
    }
    cout << '\n';
  }

  void showVector(const vec_t& x) {
    // this function prints a vector to the screen.
    using std::cout;
    cout << std::setprecision(consts::dispPrec);
    for (auto x_i : x) {
      cout << '\n' << std::setw(consts::dispPrec+2)
           << std::right << std::showpoint << x_i;
    }
    cout << '\n';
  }

  void swapRows(mat_t& A_OUT, vec_t& b_OUT, const index r1, const index r2) {
    // This function swaps row1 and row2 of both the matrix A
    // and vector b in the system Ax=b.
    double placeholder{};
    // swap rows in matrix
    for (index k{}; k < A_OUT[0].size(); ++k) {
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
    index m{    A.size() };
    index n{ A[0].size() };
    mat_t B(n,vec_t(m));
    for (index i{}; i < n; ++i) {
      for (index j{}; j < m; --j) {
        B[i][j] = A[j][i];
      }
    }
    return B;
  }

  vec_t updateVector(const vec_t& vec1, const double scalar, const vec_t& vec2) {
    // This function returns the linear update x = vec1 + scalar*vec2.
    // Here we choose to return by value (rather than pass x_OUT by reference)
    // in order to increase readability.
    index m{vec1.size()};
    vec_t  x(m);
    for (index i{}; auto vec2_i : vec2) {
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
    index m{A.size()};
    index n{A[0].size()};
    // write into the file
    outf << m << ' ';
    outf << n << '\n';
    std::setprecision(consts::writePrec);
    for (index i{}; i < m; ++i) {
      for (index j{}; j < n; ++j) {
        outf << std::setw(consts::writePrec+2)
             << std::left << std::showpoint << A[i][j] << ' ';
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
    index m{x.size()};
    outf << m << '\n';
    // write into the file
    std::setprecision(consts::writePrec);
    for (auto x_i : x) {
      outf << std::setw(consts::writePrec+2)
           << std::left << std::showpoint << x_i << '\n';
    }
    // when outf goes out of scope, the ofstream destructor will close the file
  }

}
