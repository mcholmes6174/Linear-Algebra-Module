/* Written by Matthew Holmes
 * November 2021
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
 *
 */
#include "ConjugateGrad.h"
#include "LinAlgTools.h"
#include <cmath>
#include <iostream>

void updateVector(double x[], const double vec1[], const double scalar,
                              const double vec2[], const int m=consts::m) {
  // this function performs the linear update x = vec1 + scalar*vec2
  for (int i{}; i < m; ++i) {
    x[i] = vec1[i] + scalar*vec2[i];
  }
}

void testMaxIter(int num_iter, int n=consts::n) {
  // this function checks whether the # of iterations executed in the
  // C-G method exceeds the theoretical maximum before convergence is
  // guaranteed, and throws an error if it does. It also prints every nth
  // iteration as the program executes.
  try {
    if (num_iter > consts::m) {
      throw "No convergence in theoretical maximum # iter. Exiting program.";
    }
    else {
      if (num_iter % n == 0) std::cout << "\nIteration number\t" << num_iter;
    }
  }
  catch (const char* exception) {
    std::cerr << "\nError: " << exception << '\n' << '\n';
    std::exit(0);
  }
}

void basicCG(const double A[][consts::n], const double b[], double x[]) {
  // this function executes the basic conjugate-gradient method

  // first, check that the matrix is symmetric
  tools::checkSymm(A);

  // initialize residual r, initial conjugate vector p, and residual norm
  // we assume that the initial guess to the solution is x=0
  double r[consts::m]{};
  tools::makeVecCopy(r,b);
  double p[consts::m]{};
  tools::makeVecCopy(p,b);
  double r_norm{tools::getNorm(r)};

  // iterate until the residual falls within the desired accuracy
  int iter_counter{0};
  // the next line is to get "step" to use in testMaxIter()
  //  -> we take "step" to be the smallest power of ten 10^k s.t. 10^k <= m
  //  -> to compute this, we find that 10^k <= m iff k <= log(m), so we take
  //  -> k = floor( log(m) ) and then step = 10^k
  int step{static_cast<int>(std::pow(10.0,std::floor(std::log10(consts::m))))};
  while (r_norm > consts::tolerance) {

    // store y = A*b to use in new residual and constant beta
    double y[consts::m]{};
    tools::matVecMul(A,p,y);
    // calculate alpha to minimize obj. func. along direction p
    double alpha{ tools::innerProd(p,r)/tools::innerProd(p,y) };
    // update solution vector x
    updateVector(x,x,alpha,p);
    // update residual vector to be
    // r = r - alpha*y = r - alpha*A*p = b - A*x - alpha*A*p = b - A*(x+alpha*p)
    updateVector(r,r,-alpha,y);
    // update the residual norm
    r_norm = tools::getNorm(r);
    // calculate beta to obtain next conjugate vector p
    double beta{ -tools::innerProd(r,y)/tools::innerProd(p,y) };
    // get next conjugate vector
    updateVector(p,r,beta,p);
    // keep track of # of iterations and throw error once they exceed the
    // theoretical maximum
    testMaxIter(iter_counter,step);
    iter_counter += 1;

  }
  std::cout << "\nConvergence acheived in " << iter_counter << " iterations!\n";
}

void smartCG(const double A[][consts::n], const double b[], double x[]) {
  // this function executes the smart conjugate-gradient method

  // first, check that the matrix is symmetric
  tools::checkSymm(A);

  // initialize residual r, initial conjugate vector p, and (residual norm)^2
  // we assume that the initial guess to the solution is x=0
  double r[consts::m]{};
  tools::makeVecCopy(r,b);
  double p[consts::m]{};
  tools::makeVecCopy(p,b);
  double rTr{tools::innerProd(r,r)}; // norm^2 = < r, r >

  // iterate until the residual falls within the desired accuracy
  int iter_counter{0};
  // the next line is to get "step" to use in testMaxIter()
  //  -> we take "step" to be the smallest power of ten 10^k s.t. 10^k <= m
  //  -> to compute this, we find that 10^k <= m iff k <= log(m), so we take
  //  -> k = floor( log(m) ) and then step = 10^k
  int step{static_cast<int>(std::pow(10.0,std::floor(std::log10(consts::m))))};
  while (std::sqrt(rTr) > consts::tolerance) {

    // store y = A*b to use in new residual and constant beta
    double y[consts::m]{};
    tools::matVecMul(A,p,y);
    // calculate alpha to minimize obj. func. along direction p
    double alpha{ rTr/tools::innerProd(p,y) };
    // update solution vector x
    updateVector(x,x,alpha,p);
    // update residual vector to be
    // r = r - alpha*y = r - alpha*A*p = b - A*x - alpha*A*p = b - A*(x+alpha*p)
    updateVector(r,r,-alpha,y);
    // update (residual norm)^2
    double rTr_new{tools::innerProd(r,r)};
    // calculate beta to obtain next conjugate vector p
    double beta{ rTr_new/rTr };
    // get next conjugate vector
    updateVector(p,r,beta,p);
    // update the old inner product with the new
    rTr = rTr_new;
    // keep track of # of iterations and throw error once they exceed the
    // theoretical maximum
    testMaxIter(iter_counter,step);
    iter_counter += 1;

  }
  std::cout << "\nConvergence acheived in " << iter_counter << " iterations!\n";
}

void smartPreCondCG(const double A[][consts::n], const double b[], double x[]) {
  // this function executes the smart conjugate-gradient
  // method with pre-conditioning

  // first, check that the matrix is symmetric
  tools::checkSymm(A);

  // we assume that the initial guess to the solution is x=0
  // initialize residual r
  double r[consts::m]{};
  tools::makeVecCopy(r,b);
  // compute z using preconditoner matrix
  double z[consts::m]{};
  tools::diagPreCond(A,r,z);
  // initialize conjugate vector p
  double p[consts::m]{};
  tools::makeVecCopy(p,z);
  // initialize "error" E = rTz = < r, z >
  double E{tools::innerProd(r,z)};

  // iterate until the residual falls within the desired accuracy
  int iter_counter{0};
  // the next line is to get "step" to use in testMaxIter()
  //  -> we take "step" to be the smallest power of ten 10^k s.t. 10^k <= m
  //  -> to compute this, we find that 10^k <= m iff k <= log(m), so we take
  //  -> k = floor( log(m) ) and then step = 10^k
  int step{static_cast<int>(std::pow(10.0,std::floor(std::log10(consts::m))))};
  while (std::sqrt(E) > consts::tolerance) {

    // store y = A*b to use in new residual and constant beta
    double y[consts::m]{};
    tools::matVecMul(A,p,y);
    // calculate alpha to minimize obj. func. along direction p
    double alpha{ E/tools::innerProd(p,y) };
    // update solution vector x
    updateVector(x,x,alpha,p);
    // update residual vector to be
    // r = r - alpha*y = r - alpha*A*p = b - A*x - alpha*A*p = b - A*(x+alpha*p)
    updateVector(r,r,-alpha,y);
    // update z vector
    tools::diagPreCond(A,r,z);
    // update (residual norm)^2
    double E_new{tools::innerProd(r,z)};
    // calculate beta to obtain next conjugate vector p
    double beta{ E_new/E };
    // get next conjugate vector
    updateVector(p,z,beta,p);
    // update the old inner product with the new
    E = E_new;
    // keep track of # of iterations and throw error once they exceed the
    // theoretical maximum
    testMaxIter(iter_counter,step);
    iter_counter += 1;

  }
  std::cout << "\nConvergence acheived in " << iter_counter << " iterations!\n";
}
