/* Written by Matthew Holmes
 * December 2021
 *
 * This program is the driver for my first significant C++ code. The driver
 * calls routines defined in the following files:
 *
 * ConjugateGrad.cpp
 * GuassianElim.cpp
 * GMRES.cpp
 * InOut.cpp
 * LinAlgToolkit.cpp
 * LSQ.cpp
 *
 * in order to solve a linear system Ax=b using a method chosen by the user
 * at runtime. The matrix A and vector b that define the system are to be
 * read in from two separate files (e.g., mat_A.dat and vec_b.dat). The names
 * of these two files may be provided by the user as command line arguments, but
 * if omited, the program will then prompt the user to enter the filenames.
 *
 * After reading in the contents of the two files, the program prompts the user
 * to choose a numerical method to apply. Once a selection is made, the program
 * shows the results and writes the solution vector to a file named x_soln.dat.
 *
 * The Matrix.h and Vector.h files contain the declarations of the Matrix and
 * Vector user-defined classes respectively, while the corresponding .cpp files
 * contains their member function implementations.
 *
 */
#include "ConjugateGrad.h" // Conjugate-Gradient routines
#include "GaussianElim.h"  // Gaussian Elimination routines
#include "GMRES.h"         // GMRES routines
#include "InOut.h"         // io::  routines in user-defined namespace
#include "LinAlgToolkit.h" // tlk:: routines in user-defined namespace
#include "LSQ.h"           // QR decomp and back-sub routines for least squares
#include "Matrix.h"        // user-defined Matrix class
#include "Types.h"         // type aliases for size_t, mat_t and vec_t
#include "Vector.h"        // user-defined Vector class
#include <iostream>        // std::cout and std::cin
#include <string>          // std::string

int main(int argc, char* argv[]) { // command line args will be 2 filenames
  using std::cout;
  cout << "*******************************************************************";
  cout << "\nWelcome to the Linear Algebra Module.";
  cout << "\nHere we will solve a linear system Ax=b.\n";

  // to create files containing a symmetric system for future use in CG method
  // inout::generateSymmSystem(100,100);

  // to get filenames, we first check for command line arguments
  std::string filename_A{};
  std::string filename_b{};
  std::string filename_x{"x_soln.dat"};
  if (argc < 2) {
    // if not provided, we ask user to specify filename of matrix A
    // and vector b to be used in system
    filename_A = inout::askFileChoice('A');
    filename_b = inout::askFileChoice('b');
  }
  else {
    // otherwise, use given command line arguments
    filename_A = argv[1];
    filename_b = argv[2];
  }

  // read in matrix size from file
  index m{ tlk::readNthVal(filename_A,0) };
  index n{ tlk::readNthVal(filename_A,1) };

  // create Matrix object and load with values from file
  Matrix A{m,n};
  A.load(filename_A);

  // create Vector object and load with values from file
  Vector b{m};
  b.load(filename_b);

  // create Vector object and leave it zero-initialized for future use
  Vector x{n};

  // report the size of the matrix to the user
  cout << "\nThe matrix contained in file " << filename_A;
  cout << " has size " << m << " by " << n << '\n';

  // show the system to the user
  cout << "\nHere is the matrix A:" << A;
  cout << "\nHere is the vector b:" << b;

  // ask user to choose numerical method
  int user_input{ inout::askMethodChoice(m,n) };

  // define POD enumerator type to make the following switch cases more readable
  enum Method {
    method_GaussElim,
    method_basicCG,
    method_smartCG,
    method_smartPreCondCG,
    method_GMRES,
    method_HouseholderLSQ
  };

  // apply the chosen method
  switch(user_input) {
    case method_GaussElim: {

      pivotElim(A,b);
      backSub(A,b,x);
      cout << "\nMatrix A after Gaussian Elimination with pivoting:" << A;
      cout << "\nVector b after Gaussian Elimination with pivoting:" << b;
      break;

    }
    case method_basicCG: {

      basicCG(A,b,x);
      break;

    }
    case method_smartCG: {

      smartCG(A,b,x);
      break;

    }
    case method_smartPreCondCG: {

      smartPreCondCG(A,b,x);
      break;

    }
    case method_GMRES: {

      applyGMRES(A,b,x);
      break;

    }
    case method_HouseholderLSQ: {

      Vector v_diag{m};
      decompQR(A,b,v_diag);
      solveTriLSQ(A,b,x);

    }
  }

  // show the solution vector
  cout << "\nThe solution vector x to the system is:" << x;

  // write solution vector to file
  x.write(filename_x);

  // re-read system from file in order to compute error
  A.load(filename_A);
  b.load(filename_b);

  // compute and show the error Ax-b to verify solution
  cout << "\nThe Euclidean norm of the error vector err = Ax-b is given by:";
  cout << '\n' << ((A*x)-b).getNorm() << '\n';

  cout << "\nGoodbye\n";
  cout << "*******************************************************************";
  cout << '\n';
  return 0;
}
