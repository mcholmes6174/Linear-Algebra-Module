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
#include <cmath>       // for std::pow() and std::sqrt() and std::abs()
#include <fstream>     // to do file I/O
#include <iomanip>     // to override the default precision shown by std::cout
#include <iostream>
#include <string>

namespace tools {

  int dispPrec{8};   // precision with which to display double values on screen
  int writePrec{16}; // precision with which to write double values to .dat file

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

  void checkSymm(const double A[][consts::n], const int m) {
    // this function checks whether a given square matrix is symmetric, and
    // throws an error if it is not symmetric.
    for (int i{}; i < m; ++i) {
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

  void diagPreCond(const double M[][consts::n], const double r[], double z[],
                   const int m) {
    // this function performs diagonal preconditioning for C-G method
    // via the multiplication z = M^(-1)*r, where M=diag(A) is the diagonal
    // preconditoner matrix. The original matrix A should be passed as the first
    // argument.
    for (int i{}; i < m; ++i) {
      z[i] = (1.0/M[i][i])*r[i];
    }
  }

  void getError(const double A[][consts::n], const double   x[],
                const double b[],                  double err[],
                const int m, const int n) {
    // this function computes the error in the solution x to the system Ax=b
    // once a solution has been found. It simply returns err = Ax-b.
    // need to forward declare matVecMul here
    void matVecMul(const double A[][consts::n], const double vecIN[],
                         double vecOUT[], const int, const int);
    matVecMul(A, x, err, m, n); // multiply A*x
    for (int i{}; i < m; ++i) {
      err[i] -= b[i];             // subtract b
    }
  }

  double getNorm(const double x[], const int m) {
    // this function returns the Euclidean norm of a vector x.
    double norm{};
    for (int i{}; i < m; ++i) {
      norm += std::pow(x[i],2.0);
    }
    return std::sqrt(norm);
  }

  double innerProd(const double x[], const double y[], const int m) {
    // this function returns the inner product of two vectors x and y.
    double result{};
    for (int i{}; i < m; ++i) {
      result += x[i]*y[i];
    }
    return result;
  }

  void matMul(const double A[][consts::n], const double B[][consts::n],
                      double C[][consts::n], const int m, const int n) {
    // this function performs the matrix multiplication C=AB.
    for (int i{}; i < m; ++i) {
      for (int j{}; j < n; ++j) {
        for (int k{}; k < n; ++k) {
          C[i][j] += A[i][k] * B[k][j];
        }
      }
    }
  }

  void makeVecCopy(double x_copy[], const double x[], const int m) {
    // this function copies the contents of vector x into vector x_copy.
    for (int i{}; i < m; ++i) {
      x_copy[i] = x[i];
    }
  }

  void matVecMul(const double A[][consts::n], const double x[], double y[],
                   const int m, const int n) {
    // this function performs the matrix-vector multiplication y=Ax.
    for (int i{}; i < m; ++i) {
      for (int j{}; j < n; ++j) {
        y[i] += A[i][j] * x[j];
      }
    }
  }

  void readMatrix(const std::string filename, double A[][consts::n],
                  const int m, const int n) {
    // this function reads the contents of "filename" and stores the data in
    // the matrix A.
    // we'll read from a file called "filename"
    std::ifstream inf{filename};
    // in case we cannot open the input file stream
    if (!inf) {
      std::cerr << "\nError: " << filename
                               << " could not be opened for reading\n\n";
      return;
    }
    // store the first m*n pieces of data in the matrix A
    for (int k{}; k < m*n; ++k) {
      int i{ static_cast<int>(std::floor(k/m)) };
      int j{ k % n };
      std::string strVal;
      inf >> strVal;
      A[i][j] = std::stod(strVal);
    }
    if (inf) {
      std::cerr << "\nWarning: Contents of " << filename
                                             << " not read in entirety.\n\n";
    }
  }

  void showMatrix(const double A[][consts::n], const int m, const int n) {
    // this function prints a matrix to the screen.
    using std::cout;
    cout << std::setprecision(dispPrec);
    for (int i{}; i < m; ++i) {
      cout << '\n';
      for (int j{}; j < n; ++j) {
        cout << std::setw(dispPrec+2) << std::right
                                      << std::showpoint << A[i][j] << '\t';
      }
    }
    cout << '\n';
  }

  void showVector(const double x[], const int m) {
    // this function prints a vector to the screen.
    using std::cout;
    cout << std::setprecision(dispPrec);
    for (int i{}; i < m; ++i) {
      cout << '\n' << std::setw(dispPrec+2) << std::right
                                            << std::showpoint << x[i];
    }
    cout << '\n';
  }

  void swapRows(double A[][consts::n], double b[],
                const int row1, const int row2, const int n) {
    // This function swaps row1 and row2 of both the matrix A and vector b in
    // the system Ax=b
    double placeholder{};
    // swap rows in matrix
    for (int k{}; k < n; ++k) {
      placeholder  = A[row1][k];
      A[row1][k]   = A[row2][k];
      A[row2][k]   = placeholder;
    }
    // swap rows in vector
    placeholder = b[row1];
    b[row1]     = b[row2];
    b[row2]     = placeholder;
  }

  void writeVector(const std::string filename, const double x[], const int m) {
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
    for (int i{}; i < m; ++i) {
      outf << std::setw(writePrec+2) << std::right
                                     << std::showpoint << x[i] << '\n';
    }
    // when outf goes out of scope, the ofstream destructor will close the file
  }

}
