/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains routines used in the initial I/O. This module was
 * created in order to remove clutter from the driver routine.
 *
 */
 #include "InOut.h"
 #include "LinAlgToolkit.h"
 #include "Types.h"
 #include <fstream>
 #include <iostream>
 #include <string>
 #include <vector>

namespace inout {

  using std::size_t;

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
      else cout << "\nNot a valid input, provide an integer value in [0," << n
                << ").\n";
    }
    return user_input;
  }

  void generateSymmSystem(const size_t m, const size_t n) {
    mat_t A(m,vec_t(n));
    vec_t b(n);
    double fill{1.0};
    // fill the arrays with nonzero values
    for (size_t i{}; i < m; ++i) {
      for (size_t j{}; j < n; ++j) {
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
    tlk::writeMatrix("mat_symm.dat",A);
    tlk::writeVector("vec_symm.dat",b);
  }

  size_t askSubspaceDim(const size_t m) {
    // This function is called when the GMRES method is selected. It prompts the
    // user to enter the dimension of the Krylov subspace that the solution
    // should be minimized over.
    size_t user_input{};
    bool   not_valid{true};
    while (not_valid) {
      using std::cout;
      cout << "\nPlease choose dimension for Krylov subspace:\n";
      std::cin >> user_input;

      // make sure std::cin is robust
      if (std::cin.fail()) { // to handle a failed extraction
        std::cin.clear();    // failure causes user_input to be zero-initialized
        user_input = m+1;    // so we set to m+1 to cause not_valid to remain true
      }
      ignoreLine();          // to remove any extraneous input

      // check to see whether input is a valid integer
      if (user_input <= m) not_valid = false;
      else cout << "\nNot a valid input, provide an integer value in [0," << m
                << "].\n";
    }
    return user_input;
  }

}
