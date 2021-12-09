/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains all the routines necessary to perform the Generalized
 * Minimal Residual (GMRES) method for solving a square linear system Ax=b.
 *
 * NOTE: GMRES can be applied to ANY square system provided A is invertible.
 *
 */
 #include "Constants.h"
 #include "GMRES.h"
 #include "LinAlgToolkit.h"
 #include "LSQ.h"
 #include "Matrix.h"
 #include "Types.h"
 #include "Vector.h"
 #include <iostream>
 #include <string>

void applyArnoldi(const Matrix& A, Matrix& Q_OUT, Matrix& H_OUT, const index j) {
  // This function applies the jth step of the Arnoldi Algorithm in order to
  // compute an orthonormal basis for the Krylov subspace K_n(A,b), where the
  // vector b is implicitly contained in the first column of Q as
  // q_1 = r_0/beta where r_0 = b - A*x_0 and beta = ||r_0||_2.
  //
  // The matrix Q is modified directly, with its columns containing the
  // orthonormal basis on exit. The matrix H is modified directly as well,
  // becoming the n+1 by n Hessenberg version of A on exit.

  // get sizes
  index m{ A.size(0) };

  // declare and zero-initialize column vectors
  Vector q_next{m};                    // q_{j+1}
  Vector q_j{m};                       // q_j
  tlk::nthColumn("grab",Q_OUT,q_j,j);  // get q_j from Q_OUT
  tlk::matVecMul(A,q_j,q_next);        // q_{j+1} = A*q_j

  // loop under index j and modify q to ensure mutual orthogonality
  for (index i{}; i <= j; ++i) {
    H_OUT.set(i,j) = tlk::innerProd(q_next,q_j);  // < q_{j+1}, q_i >
    Vector q_i{m};
    tlk::nthColumn("grab",Q_OUT,q_i,i);
    q_next = tlk::updateVector(q_next,-H_OUT.get(i,j),q_i); // q_{j+1} += -h_{i,j}*q_i
  }

  // normalize the column vector q
  H_OUT.set(j+1,j) = q_next.getNorm();             // get norm of q_{i+j}
  if (std::abs(H_OUT.get(j+1,j)) < consts::epsilon) {  // if zero, stop
    std::cerr << "Divide by zero warning.";
    // std::exit(0);
  }
  q_next = tlk::makeVecCopy(q_next, 1.0/H_OUT.get(j+1,j) );

  /// store q_{j+1} in (j+1)th column of Q
  tlk::nthColumn("give",Q_OUT,q_next,j+1);

}

void multiplyByQ(const Matrix& Q, const Vector& y, Vector& x_OUT) {
  // this function performs the multiplication x = Q*y

  // get size of m x n system x_n = Q_n*y_n
  index m{ Q.size(0)   };
  index n{ Q.size(1)-1 };

  // multiply
  for (index i{}; i < m; ++i) {
    for (index j{}; j < n; ++j) {
      x_OUT.set(i) += Q.get(i,j) * y.get(j);
    }
  }

}

void applyGMRES(Matrix& A_OUT, Vector& b_OUT, Vector& x_OUT) {
  // This function executes the GMRES method to solve the square linear system
  // Ax=b.

  // first, we get the size of the system
  const index m{A_OUT.size(0)};
  const index max_dim{A_OUT.size(1)};

  // and we begin the iterative process with
  index n{1};
  // which will specify the dimension of the Krylov subspace we minimize over

  // We declare and initialize vector r_0, scalar beta, and matries Q_n and H_n
  Vector r{m};
  r =    tlk::makeVecCopy(b_OUT);
  double beta{ r.getNorm() };
  Matrix Q_np1{ m,  n+1};
  Matrix H_n  {n+1,  n };

  // store (1/beta)*r in first column of Q_n (a hassle b/c C is row major)
  // this is the first vector in our orthonormal basis
  for (index i{}; i < r.size(); ++i) {
    Q_np1.set(i,0) = (1.0/beta)*r.get(i);
  }

  // now, we iterate over n = 1, 2, 3, ..., A.size() until the desired tolerance
  // is acheived or we reach n = A[0].size()
  while (beta > consts::epsilon) {

    // for debugging
    std::cout << "\nn = " << n << '\n';
    // std::cout << "\nQ_{n+1}\n";
    // tlk::showMatrix(Q_np1);
    // std::cout << "\nH\n";
    // tlk::showMatrix(H_n);

    // apply nth step of Arnoldi Algorithm
    applyArnoldi(A_OUT, Q_np1, H_n, n-1);

    // for debugging
    // std::cout << "\nQ_{n+1}\n";
    // tlk::showMatrix(Q_np1);
    // std::cout << "\nH\n";
    // tlk::showMatrix(H_n);

    // find y that minimizes ||beta*e1 - H_n*y||Vector y(n);
    Vector v_diag{n};
    Vector beta_e1{n+1};
    beta_e1.set(0) = beta;
    decompQR(H_n, beta_e1, v_diag);

    // for debugging
    // std::cout << "\nv_diag\n";
    // tlk::showVector(v_diag);
    // std::cout << "\nbeta*e1\n";
    // tlk::showVector(beta_e1);
    // std::cout << "\nH\n";
    // tlk::showMatrix(H_n);

    Vector y{n};
    solveTriLSQ(H_n, beta_e1, y);

    // for debugging
    // std::cout << "\ny\n";
    // tlk::showVector(y);

    // take x_n = x_0 + Q_n*y = Q_n*y
    multiplyByQ(Q_np1,y,x_OUT);

    // for debugging
    std::cout << "\nx_n\n";
    x_OUT.show();

    // check whether we have reached the theoretical max
    if (n < max_dim) {
      // if not, increment n and continue
      n += 1;
    }
    else {
      std::cerr << "\nMaximum number of iterations reached w/o convergence.";
      std::cerr << "\nExiting the iterative process.\n";
      break;
    }

    // resize Q_{n+1} and H_n for next iteration
    Q_np1.resize(m,n+1);
    H_n.resize(n+1,n);

  }

}
