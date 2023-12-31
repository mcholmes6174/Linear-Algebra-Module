/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains all the routines necessary to perform
 *
 * 1) the basic Conjugate-Gradient
 * 2) the smart Conjugate-Gradient
 * 3) the smart preconditioned Conjugate-Gradient
 *
 * methods for solving a square linear system Ax=b.
 *
 * NOTE: The conjugate-gradient methods require the matrix A to
 *       be symmetric AND strictly positive-definite.
 */
#include "ConjugateGrad.h"
#include "LinAlgToolkit.h"
#include "Matrix.h"
#include "Types.h"
#include "Vector.h"
#include <cmath>
#include <iostream>
#include <vector>

void testMaxIter(int num_iter, index max_iter, int step) {
  // this function checks whether the # of iterations executed in the
  // C-G method exceeds the theoretical maximum before convergence is
  // guaranteed, and throws an error if it does. It also prints every nth
  // iteration as the program executes.
  try {
    if (num_iter > static_cast<int>(max_iter)) {
      throw "No convergence in theoretical maximum # iter. Exiting program.";
    }
    else {
      if (num_iter % step == 0) std::cout << "\nIteration number\t" << num_iter;
    }
  }
  catch (const char* exception) {
    std::cerr << "\nError: " << exception << '\n' << '\n';
    std::exit(0);
  }
}

void basicCG(const Matrix& A, const Vector& b, Vector& x_OUT) {
  // this function executes the basic conjugate-gradient method

  // first, check that the matrix is symmetric
  tlk::checkSymm(A);

  // then get size of system
  const index m{A.size(0)};

  // initialize residual r, initial conjugate vector p, and residual norm
  // we assume that the initial guess to the solution is x=0
  Vector r{m};
  r = b;
  Vector p{m};
  p = b;
  double r_norm{ r.getNorm() };

  // iterate until the residual falls within the desired accuracy
  int iter_counter{0};
  // the next line is to get "step" to use in testMaxIter()
  //  -> we take "step" to be the smallest power of ten 10^k s.t. 10^k <= m
  //  -> to compute this, we find that 10^k <= m iff k <= log(m), so we take
  //  -> k = floor( log(m) ) and then step = 10^k
  int step{static_cast<int>(std::pow(10.0,std::floor(std::log10(m))))};
  while (r_norm > consts::tolerance) {

    // store y = A*b to use in new residual and constant beta
    Vector y{m};
    y = A*p;
    // calculate alpha to minimize obj. func. along direction p
    double alpha{ tlk::innerProd(p,r)/tlk::innerProd(p,y) };
    // update solution vector x
    x_OUT = x_OUT + alpha*p;
    // update residual vector to be
    // r = r - alpha*y = r - alpha*A*p = b - A*x - alpha*A*p = b - A*(x+alpha*p)
    r = r - alpha*y;
    // update the residual norm
    r_norm = r.getNorm();
    // calculate beta to obtain next conjugate vector p
    double beta{ -tlk::innerProd(r,y)/tlk::innerProd(p,y) };
    // get next conjugate vector
    p = r + beta*p;
    // keep track of # of iterations and throw error once they exceed the
    // theoretical maximum
    testMaxIter(iter_counter,m,step);
    iter_counter += 1;

  }
  std::cout << "\nConvergence acheived in " << iter_counter << " iterations!\n";
}

void smartCG(const Matrix& A, const Vector& b, Vector& x_OUT) {
  // this function executes the smart conjugate-gradient method

  // first, check that the matrix is symmetric
  tlk::checkSymm(A);

  // then get size of system
  index m{A.size(0)};

  // initialize residual r, initial conjugate vector p, and (residual norm)^2
  // we assume that the initial guess to the solution is x=0
  Vector r{m};
  r = b;
  Vector p{m};
  p = b;
  double rTr{tlk::innerProd(r,r)}; // norm^2 = < r, r >

  // iterate until the residual falls within the desired accuracy
  int iter_counter{0};
  // the next line is to get "step" to use in testMaxIter()
  //  -> we take "step" to be the smallest power of ten 10^k s.t. 10^k <= m
  //  -> to compute this, we find that 10^k <= m iff k <= log(m), so we take
  //  -> k = floor( log(m) ) and then step = 10^k
  int step{static_cast<int>(std::pow(10.0,std::floor(std::log10(m))))};
  while (std::sqrt(rTr) > consts::tolerance) {

    // store y = A*b to use in new residual and constant beta
    Vector y{m};
    y = A*p;
    // calculate alpha to minimize obj. func. along direction p
    double alpha{ rTr/tlk::innerProd(p,y) };
    // update solution vector x
    x_OUT = x_OUT + alpha*p;
    // update residual vector to be
    // r = r - alpha*y = r - alpha*A*p = b - A*x - alpha*A*p = b - A*(x+alpha*p)
    r = r - alpha*y;
    // update (residual norm)^2
    double rTr_new{tlk::innerProd(r,r)};
    // calculate beta to obtain next conjugate vector p
    double beta{ rTr_new/rTr };
    // get next conjugate vector
    p = r + beta*p;
    // update the old inner product with the new
    rTr = rTr_new;
    // keep track of # of iterations and throw error once they exceed the
    // theoretical maximum
    testMaxIter(iter_counter,m,step);
    iter_counter += 1;

  }
  std::cout << "\nConvergence acheived in " << iter_counter << " iterations!\n";
}

void smartPreCondCG(const Matrix& A, const Vector& b, Vector& x_OUT) {
  // this function executes the smart conjugate-gradient
  // method with pre-conditioning

  // first, check that the matrix is symmetric
  tlk::checkSymm(A);

  // then get size of system
  index m{A.size(0)};

  // we assume that the initial guess to the solution is x=0
  // initialize residual r
  Vector r{m};
  r = b;
  // compute z using preconditoner matrix
  Vector z{m};
  tlk::diagPreCond(A,r,z);
  // initialize conjugate vector p
  Vector p{m};
  p = z;
  // initialize "error" E = rTz = < r, z >
  double E{tlk::innerProd(r,z)};

  // iterate until the residual falls within the desired accuracy
  int iter_counter{0};
  // the next line is to get "step" to use in testMaxIter()
  //  -> we take "step" to be the smallest power of ten 10^k s.t. 10^k <= m
  //  -> to compute this, we find that 10^k <= m iff k <= log(m), so we take
  //  -> k = floor( log(m) ) and then step = 10^k
  int step{static_cast<int>(std::pow(10.0,std::floor(std::log10(m))))};
  while (std::sqrt(E) > consts::tolerance) {

    // store y = A*b to use in new residual and constant beta
    Vector y{m};
    y = A*p;
    // calculate alpha to minimize obj. func. along direction p
    double alpha{ E/tlk::innerProd(p,y) };
    // update solution vector x
    x_OUT = x_OUT + alpha*p;
    // update residual vector to be
    // r = r - alpha*y = r - alpha*A*p = b - A*x - alpha*A*p = b - A*(x+alpha*p)
    r = r - alpha*y;
    // update z vector
    tlk::diagPreCond(A,r,z);
    // update (residual norm)^2
    double E_new{tlk::innerProd(r,z)};
    // calculate beta to obtain next conjugate vector p
    double beta{ E_new/E };
    // get next conjugate vector
    p = z + beta*p;
    // update the old inner product with the new
    E = E_new;
    // keep track of # of iterations and throw error once they exceed the
    // theoretical maximum
    testMaxIter(iter_counter,m,step);
    iter_counter += 1;

  }
  std::cout << "\nConvergence acheived in " << iter_counter << " iterations!\n";
}
