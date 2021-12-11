/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains routines used in the initial I/O. This module was
 * created in order to remove clutter from the driver routine.
 *
 */
 #include "InOut.h"
 #include "LinAlgToolkit.h"
 #include "Matrix.h"
 #include "Types.h"
 #include "Vector.h"
 #include <fstream>
 #include <iostream>
 #include <string>

namespace inout {

  int listMethods() {
    // this function lists each of the methods we have available in the module,
    // and returns "n" which is the total number of methods
    using std::cout;
    constexpr int n{6};
    std::string method_list[n]{
      "Gaussian Elimination with Pivoting",
      "Basic Conjugate-Gradient",
      "Smart Conjugate-Gradient",
      "Smart Diagonally-Preconditioned Conjugate-Gradient",
      "Generalized Minimal Residual (GMRES)",
      "Least Squares via QR decomposition (w/ Householder transformations)"
    };
    cout << '\n';
    for (int i{}; auto method : method_list) {
      cout << "\n " << i << ") " << method << '\n';
      ++i;
    }
    cout << '\n';
    return n;
  }

  void ignoreLine() {
    // this function tells std::cin to ignore all buffered characters through
    // and including the next newline character
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
  }

  std::string askFileChoice(const char A_or_b) {
    using std::cout;
    switch (A_or_b) {
      case 'A': {
        cout << "\nPlease enter name of file containing system matrix A:\n\n";
        break;
      }
      case 'b': {
        cout << "\nPlease enter name of file containing system vector b:\n\n";
        break;
      }
    }
    std::string filename{};
    bool not_valid{true};
    while (not_valid) {
      std::cin >> filename;
      // make sure std::cin is robust
      if (std::cin.fail()) { // to handle a failed extraction
        std::cin.clear();    // failure causes user_input to be zero-initialized
      }
      ignoreLine();          // to remove any extraneous input

      // check to see whether input corresponds to existing file
      std::ifstream inf{filename};
      if (!inf) cout << "\nFile not found, please try again.\n";
      else not_valid = false;
    }
    return filename;
  }

  int askMethodChoice(const index m, const index n) {
    int  user_input{};
    bool not_valid{true};
    while (not_valid) {
      using std::cout;
      cout << "\nPlease choose a method to apply:";
      int n_methods{};
      n_methods = listMethods(); // returns total # of methods
      std::cin >> user_input;

      // make sure std::cin is robust
      if (std::cin.fail()) { // to handle a failed extraction
        std::cin.clear();    // failure causes user_input to be zero-initialized
        user_input = -1;     // so we set to -1 to cause not_valid to remain true
      }
      ignoreLine();          // to remove any extraneous input

      // check to see whether input is a valid integer
      if (0 <= user_input && user_input < n_methods) not_valid = false;
      else cout << "\nNot a valid input, provide an integer value in [0," << n
                << ").\n";

      // check to see if method is compatible with system
      const bool square{m==n};
      const bool overdetermined{m>n};
      if (square && user_input == n_methods-1) {
        not_valid = true;
        cout << "\nSystem is square, please choose method accordingly.\n";
      }
      else if (overdetermined && user_input != n_methods-1) {
        not_valid = true;
        cout << "\nSystem is overdetermined,"
             << " please choose method accordingly.\n";
      }
      /************************************************************************/
      /*********** TEMPORARY CHECK WHILE GMRES IS BEING IMPLEMENTED ***********/
      /************************************************************************/
      // we comment this section out while we debug GMRES
      if (user_input == 4) {
        not_valid = true;
        cout << "\nThe GMRES method is still being developed, please choose a"
             << " different method.\n";
      }
      /************************************************************************/
      /************************************************************************/
      /************************************************************************/
    }
    return user_input;
  }

  void generateSymmSystem(const index m, const index n) {
    Matrix A{m,n};
    Vector b{m};
    double fill{1.0};
    // fill the arrays with nonzero values
    for (index i{}; i < m; ++i) {
      for (index j{}; j < n; ++j) {
        // to handle the main diagonal differently than off-diagonal
        switch (i==j) {
          case true: {
            // we create an ascending diagonal: diag(A)=[1,2,3,4,...,m]
            A.set(i,j) = i+1;
            b.set(i)   = i+1;
            break;
          }
          default: {
            // we set the off-diagonal to be constant: off_diag(A)=fill
            A.set(i,j) = fill;
            break;
          }
        }
      }
    }
    A.write("mat_symm.dat");
    b.write("mat_symm.dat");
  }

}
