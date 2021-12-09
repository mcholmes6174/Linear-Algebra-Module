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
 */
#include "ConjugateGrad.h" // Conjugate-Gradient routines
#include "GaussianElim.h"  // Gaussian Elimination routines
#include "GMRES.h"         // GMRES routines
#include "InOut.h"         // io::  routines in user-defined namespace
#include "LinAlgToolkit.h" // tlk:: routines in user-defined namespace
#include "LSQ.h"           // QR decomp and back-sub routines for least squares
#include "Matrix.h"        // Matrix user-defined class
#include "Types.h"         // type aliases mat_t and vec_t for std::vector
#include "Vector.h"        // Vector user-defined class
#include <iostream>        // std::cout and std::cin
#include <string>          // std::string
#include <vector>          // std::vector (dynamic array functionality)

void exploreOOP() {

  // here we declare a variable of class Vector called x and set its length
  index  m{5};
  Vector x{m}; // we can call x{} or x{m} b/c of our default constructor
  Vector y{};
  y.setSize(m);

  // now we initialize the vector with nonzero values
  for (index i{}; i < m; ++i) {
    x.setVal(i,i+1);
    y.setVal(i,2);
  }

  // we access the member function called show()
  x.show();
  y.show();

  // take the inner product of x and y
  double product{ tlk::innerProd2(x,y) };
  std::cout << "\nThe inner product of x and y is: " << product << '\n';

  // here we declare a variable of class Matrix called A and set its size
  Matrix A;
  index n{5};
  A.setSize(m,n);

  // here we access the member function called show()
  A.show();

  // now we initialize the vector with nonzero values
  for (index i{}; i < m; ++i) {
    for (index j{}; j < n; ++j) {
      A.setVal(i,j,i+j);
    }
  }

  // we show again
  A.show();

}

int main(int argc, char* argv[]) { // command line args will be 2 filenames
  using std::cout;
  cout << "*******************************************************************";
  cout << "\nWelcome to the Linear Algebra Module.";
  cout << "\nHere we will solve a linear system Ax=b.\n";

  // here we will experiment with object oriented programming
  exploreOOP();
  std::exit(0);

  // to create files containing a symmetric system for future use in CG method
  // inout::generateSymmSystem(100,100);

  // to get filenames, we first check for command line arguments
  std::string filename_A{};
  std::string filename_b{};
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

  // dynamically allocate m by n array and initialize using readMatrix()
  mat_t A(m,vec_t(n));
  tlk::readMatrix(filename_A,A);

  // dynamically allocate vector of length m and initialize using readVector()
  vec_t b(m);
  tlk::readVector(filename_b,b);

  // dynamically allocate and zero-initialize the solution vector x
  vec_t x(n);

  // report the size of the matrix to the user
  cout << "\nThe matrix contained in file " << filename_A
       << " has size " << m << " by " << n << '\n';

  // show the original system to the user if small enough
  bool small_system{std::max(m,n) <= 8};
  if  (small_system) {
    cout << "\nHere is the matrix A:";
    tlk::showMatrix(A);
    cout << "\nHere is the vector b:";
    tlk::showVector(b);
  }

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

      applyGMRES(A,b,x);
      break;

    }
    case method_HouseholderLSQ: {

      vec_t v_diag(m);
      decompQR(A,b,v_diag);
      solveTriLSQ(A,b,x);

    }
  }

  // show the solution vector if system is small enough
  if (small_system) {
    cout << "\nThe solution vector x to the system is:";
    tlk::showVector(x);
  }

  // write solution vector to file
  tlk::writeVector("x_soln.dat",x);

  // re-read systemfrom file in order to compute error
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
