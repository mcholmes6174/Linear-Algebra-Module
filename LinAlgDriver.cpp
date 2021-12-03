/* Written by Matthew Holmes
 * December 2021
 *
 * This program is the driver for my first significant C++ code. The driver
 * calls routines defined in the following files:
 *
 * ConjugateGrad.cpp
 * GuassianElim.cpp
 * GMRES.cpp
 * Inout.cpp
 * LinAlgToolkit.cpp
 *
 * in order to solve a square linear system Ax=b using a method chosen by the
 * user at runtime. The matrix A and vector b that define the system are to be
 * read in from two separate files (e.g., mat_A.dat and vec_b.dat). The names
 * of these two files may be provided by the user as command line arguments, but
 * if omited, the program will then prompt the user to enter the filenames.
 * After reading in the contents from the two files, the program prompts the user
 * to choose a numerical method to apply. Once a selection is made, the program
 * shows the results and writes the solution vector to a file named x_soln.dat.
 *
 */
#include "ConjugateGrad.h" // Conjugate-Gradient routines
#include "GaussianElim.h"  // Gaussian Elimination routines
#include "GMRES.h"         // GMRES routines
#include "LinAlgToolkit.h" // tlk:: routines in user-defined namespace
#include "InOut.h"         // io::  routines in user-defined namespace
#include "Types.h"         // type aliases mat_t and vec_t for std::vector
#include <iostream>        // std::cout and std::cin
#include <string>          // std::string
#include <vector>          // std::vector (dynamic array functionality)

int main(int argc, char* argv[]) { // command line args will be 2 filenames
  using std::cout;
  cout << "*******************************************************************";
  cout << "\nWelcome to the Linear Algebra Module.";
  cout << "\nHere we will solve a square linear system Ax=b.\n";

  // to create files containing a symmetric system for future use in CG method
  // inout::generateSymmSystem(100,100);

  // check for command line arguments
  std::string filename_A{};
  std::string filename_b{};
  if (argc < 2) {
    // if not provided, ask user to specify filename of matrix A
    // and vector b to be used in system
    filename_A = inout::askFileChoice('A');
    filename_b = inout::askFileChoice('b');
  }
  else {
    // otherwise, use command line arguments
    filename_A = argv[1];
    filename_b = argv[2];
  }

  // read in matrix size from file
  size_t m{ tlk::readNthVal(filename_A,0) };
  size_t n{ tlk::readNthVal(filename_A,1) };

  // dynamically allocate m by n array and initialize using readMatrix()
  mat_t A(m,vec_t(n));
  tlk::readMatrix(filename_A,A);
  // dynamically allocate vector of length m and initialize using readVector()
  vec_t b(m);
  tlk::readVector(filename_b,b);
  // declare and zero-initialize the solution vector x
  vec_t x(m);

  // report the size of the matrix to the user
  cout << "\nThe matrix contained in file " << filename_A
       << " has size " << m << " by " << n << '\n';

  // show the original system to the user if small enough
  bool small_system{std::max(m,n) <= 8};
  if (small_system) {
    cout << "\nHere is the matrix A:";
    tlk::showMatrix(A);
    cout << "\nHere is the vector b:";
    tlk::showVector(b);
  }

  // ask user to choose numerical method
  int user_input{ inout::askMethodChoice() };

  // define POD enumerator type to make switch cases more clear
  enum Method {
    method_GaussElim,
    method_basicCG,
    method_smartCG,
    method_smartPreCondCG,
    method_GMRES
  };

  // apply the chosen method
  switch(user_input) {
    case method_GaussElim: {
      pivotElim(A,b);
      if (small_system) {
        cout << "\nHere is matrix A after Gaussian Elimination with pivoting:";
        tlk::showMatrix(A);
        cout << "\nHere is vector b after Gaussian Elimination with pivoting:";
        tlk::showVector(b);
      }
      backSub(A,b,x);
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
      applyGMRES(A,b,x, inout::askSubspaceDim(m) );
      break;
    }
  }

  // show the solution vector if system is small enough
  if (small_system) {
    cout << "\nThe solution vector x to the system is:";
    tlk::showVector(x);
  }

  // write solution vector to file
  tlk::writeVector("x_soln.dat",x);

  // re-read system in order to compute error
  tlk::readMatrix(filename_A,A);
  tlk::readVector(filename_b,b);

  // compute and show the error Ax-b to verify solution
  vec_t err(m);
  tlk::getError(A,x,b,err);
  cout << "\nThe Euclidean norm of the error vector err = Ax-b is given by:";
  cout << '\n' << tlk::getNorm(err) << '\n';

  cout << "\nGoodbye\n";
  cout << "*******************************************************************";
  return 0;
}
