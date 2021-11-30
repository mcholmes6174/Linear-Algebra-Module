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
#include <string>          // for std::string
#include <iostream>        // for std::cout and std::cin

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

void generateSystem(double A[][consts::n], double A_fixed[][consts::n],
                    double b[],            double b_fixed[], double fill) {
  // fill the arrays with nonzero values
  for (int i{}; i < consts::m; ++i) {
    for (int j{}; j < consts::n; ++j) {
      // to handle the main diagonal differently than off-diagonal
      switch (i==j) {
        case true: {
          // we create an ascending diagonal: diag(A)=[1,2,3,4,...,m]
          A      [i][j] = i+1;
          A_fixed[i][j] = i+1;
          b      [i]    = i+1;
          b_fixed[i]    = i+1;
          break;
        }
        default: {
          // we set the off-diagonal to be constant: off_diag(A)=fill
          A      [i][j] = fill;
          A_fixed[i][j] = fill;
          break;
        }
      }
    }
  }
}

int main() {
  using std::cout;
  cout << "*******************************************************************";
  cout << "\nWelcome to the Linear Algebra Module.";
  cout << "\nHere we will solve a square linear system Ax=b.\n";

  // read in matrix from file
  // double test[consts::m][consts::n]{};
  // tools::readMatrix("matrix_A.dat",test);
  // tools::showMatrix(test);
  // std::exit(0);

  // declare fixed 2D array A and fixed 1D array b for the system Ax=b
  // Two instances of each will be created: one is to be altered during
  // computation, while the other is to later verify the result
  double A      [consts::m][consts::n]{};
  double A_fixed[consts::m][consts::n]{};
  double b      [consts::m]{};
  double b_fixed[consts::m]{};
  double fill{1};
  generateSystem(A,A_fixed,b,b_fixed,fill);

  // declare and zero-initialize the solution vector x
  double x[consts::m]{};

  // show the original system to the user if small enough
  bool small_system{consts::m <= 8};
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

  // compute and show the error Ax-b to verify solution
  double err[consts::m]{};
  tools::getError(A_fixed,x,b_fixed,err);
  double norm{};
  norm = tools::getNorm(err);
  cout << "\nThe Euclidean norm of the error vector err = Ax-b given by:";
  cout << '\n' << norm << '\n';

  cout << "\nGoodbye\n";
  cout << "*******************************************************************";
  return 0;
}
