/* Written by Matthew Holmes
 * November 2021
 *
 * This program is the driver for my first significant C++ code. The driver
 * calls routines defined in the GuassianElim.cpp, ConjugateGrad.cpp, and
 * GMRES.cpp files in order to solve a square linear system Ax=b using a method
 * chosen by the user at runtime.
 *
 * In addition to printing the solution to the screen, this program writes the
 * solution vector to a file named x_soln.dat. Routines from the LinAlgTools.cpp
 * file are called often to perform common operations and I/O procedures.
 *
 */
#include "Constants.h"     // for consts:: in user-defined namespace
#include "ConjugateGrad.h" // for Conjugate Gradient and GMRES routines
#include "GaussianElim.h"  // for Gaussian Elimination routines
#include "LinAlgTools.h"   // for tools:: in user-defined namespace
#include "Types.h"         // for type aliases mat_t and vec_t
#include <iostream>        // for std::cout and std::cin
#include <string>          // for std::string
#include <vector>          // for std::vector (dynamic array functionality)

int listMethods() {
  // this function lists each of the methods we have available in the module,
  // and returns "n" which is the total number of methods
  using std::cout;
  constexpr int n{5};
  std::string method_list[n]{
    "Gaussian Elimination with Pivoting and Back-Substitution",
    "Conjugate-Gradient",
    "Smart Conjugate-Gradient",
    "Smart Diagonally-Preconditioned Conjugate-Gradient",
    "Generalized Minimal Residual (GMRES)"
  };
  for (int i{}; i < n; ++i) {
    if (i == 0) cout << '\n';
    cout << "\n " << i << ") " << method_list[i] << '\n';
    if (i == n-1) cout << '\n';
  }
  return n;
}

void ignoreLine() {
  // this function tells std::cin to ignore all buffered characters through and
  // including the next newline character
  std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
}

int askMethodChoice() {
  int  user_input{};
  bool not_valid{true};
  while (not_valid) {
    using std::cout;
    cout << "\nPlease choose a method to apply:";
    int n{};
    n = listMethods(); // returns total # of methods
    std::cin >> user_input;

    // make sure std::cin is robust
    if (std::cin.fail()) { // to handle a failed extraction
      std::cin.clear();    // failure causes user_input to be zero-initialized
      user_input = -1;     // so we set to -1 to cause not_valid to remain true
    }
    ignoreLine();          // to remove any extraneous input

    // check to see whether input is a valid integer
    if (0 <= user_input && user_input < n) not_valid = false;
    else cout << "\nNot a valid input, provide an integer value in [0,5).\n";
  }
  return user_input;
}

void generateSystem(int m, int n) {
  mat_t A(m,vec_t(n));
  vec_t b(n);
  double fill{1.0};
  // fill the arrays with nonzero values
  for (int i{}; i < m; ++i) {
    for (int j{}; j < n; ++j) {
      // to handle the main diagonal differently than off-diagonal
      switch (i==j) {
        case true: {
          // we create an ascending diagonal: diag(A)=[1,2,3,4,...,m]
          A      [i][j] = i+1;
          b      [i]    = i+1;
          break;
        }
        default: {
          // we set the off-diagonal to be constant: off_diag(A)=fill
          A      [i][j] = fill;
          break;
        }
      }
    }
  }
  tools::writeMatrix("matrix_A_CG.dat",A);
  tools::writeVector("vector_b_CG.dat",b);
}

int main() {
  using std::cout;
  cout << "*******************************************************************";
  cout << "\nWelcome to the Linear Algebra Module.";
  cout << "\nHere we will solve a square linear system Ax=b.\n";

  // create file with CG system
  generateSystem(10,10);

  // set filename of matrix A and vector b to be used in system
  std::string filename_A{"matrix_A.dat"};
  std::string filename_b{"vector_b.dat"};

  // read in matrix size from file
  int m{ tools::readNthVal(filename_A,0) };
  int n{ tools::readNthVal(filename_A,1) };

  // report the size of the matrix to the user
  cout << "\nThe matrix contained in file " << filename_A
       << " has size " << m << " by " << n << '\n';

  // dynamically allocate m by n array and initialize using readMatrix()
  mat_t A(m,vec_t(n));
  A = tools::readMatrix(filename_A);
  // dynamically allocate vector of length m and initialize using readVector()
  vec_t b(m);
  b = tools::readVector(filename_b);
  // declare and zero-initialize the solution vector x
  vec_t x(m);

  // show the original system to the user if small enough
  bool small_system{std::max(m,n) <= 8};
  if (small_system) {
    cout << "\nHere is the matrix A:";
    tools::showMatrix(A);
    cout << "\nHere is the vector b:";
    tools::showVector(b);
  }

  // ask user to choose numerical method
  int user_input{ askMethodChoice() };

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
        tools::showMatrix(A);
        cout << "\nHere is vector b after Gaussian Elimination with pivoting:";
        tools::showVector(b);
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
    tools::showVector(x);
  }

  // write solution vector to file
  tools::writeVector("x_soln.dat",x);

  // re-read system in order to compute error
  A = tools::readMatrix(filename_A);
  b = tools::readVector(filename_b);

  // compute and show the error Ax-b to verify solution
  vec_t err(m);
  err = tools::getError(A,x,b);
  cout << "\nThe Euclidean norm of the error vector err = Ax-b is given by:";
  cout << '\n' << tools::getNorm(err) << '\n';

  cout << "\nGoodbye\n";
  cout << "*******************************************************************";
  return 0;
}
