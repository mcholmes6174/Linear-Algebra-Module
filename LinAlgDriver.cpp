/* Written by Matthew Holmes
 * November 2021
 *
 * This program is the driver for my first significant C++ code. The driver
 * calls routines defined in the GuassianElim.cpp, ConjugateGrad.cpp, Inout.cpp,
 * and GMRES.cpp files in order to solve a square linear system Ax=b using a
 * method chosen by the user at runtime.
 *
 * In addition to printing the solution to the screen, this program writes the
 * solution vector to a file named x_soln.dat. Routines from the LinAlgTools.cpp
 * file are called often to perform common operations and I/O procedures.
 *
 */
#include "ConjugateGrad.h" // for Conjugate Gradient and GMRES routines
#include "GaussianElim.h"  // for Gaussian Elimination routines
#include "LinAlgToolkit.h" // for tlk:: in user-defined namespace
#include "InOut.h"         // for io::  in user-defined namespace
#include "Types.h"         // for type aliases mat_t and vec_t
#include <iostream>        // for std::cout and std::cin
#include <string>          // for std::string
#include <vector>          // for std::vector (dynamic array functionality)

int main() {
  using std::cout;
  cout << "*******************************************************************";
  cout << "\nWelcome to the Linear Algebra Module.";
  cout << "\nHere we will solve a square linear system Ax=b.\n";

  // to create files containing a symmetric system for us in CG method
  // inout::generateSymmSystem(100,100);

  // ask user to specify filename of matrix A and vector b to be used in system
  std::string filename_A{ inout::askFileChoice('A') };
  std::string filename_b{ inout::askFileChoice('b') };

  // read in matrix size from file
  size_t m{ tlk::readNthVal(filename_A,0) };
  size_t n{ tlk::readNthVal(filename_A,1) };

  // dynamically allocate m by n array and initialize using readMatrix()
  mat_t A(m,vec_t(n));
  A = tlk::readMatrix(filename_A);
  // dynamically allocate vector of length m and initialize using readVector()
  vec_t b(m);
  b = tlk::readVector(filename_b);
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
  A = tlk::readMatrix(filename_A);
  b = tlk::readVector(filename_b);

  // compute and show the error Ax-b to verify solution
  vec_t err(m);
  err = tlk::getError(A,x,b);
  cout << "\nThe Euclidean norm of the error vector err = Ax-b is given by:";
  cout << '\n' << tlk::getNorm(err) << '\n';

  cout << "\nGoodbye\n";
  cout << "*******************************************************************";
  return 0;
}
