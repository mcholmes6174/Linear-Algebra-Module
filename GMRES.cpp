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
 #include <cassert>
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
  // becoming an n+1 by n Hessenberg matrix satisfying A*Q_n = Q_{n+1}*H_n.

  // get column length
  index m{ A.size(0) };

  // declare and zero-initialize column vectors
  Vector q_next{m};            // q_{j+1}
  Vector q_j{m};               // q_j
  q_j = Q_OUT.getCol(j);       // get q_j from Q
  q_next = A*q_j;              // set q_{j+1} = A*q_j

  // loop under index j and modify q to ensure mutual orthogonality
  for (index i{}; i <= j; ++i) {
    H_OUT(i,j) = tlk::innerProd(q_next,q_j);  // < q_{j+1}, q_i >
    Vector q_i{m};
    q_i = Q_OUT.getCol(i);
    q_next = q_next - H_OUT(i,j)*q_i; // set q_{j+1} -= h_{i,j}*q_i
  }

  // normalize the column vector q
  H_OUT(j+1,j) = q_next.getNorm();                 // get norm of q_{i+j}
  if (std::abs(H_OUT(j+1,j)) < consts::epsilon) {  // if zero, stop
    std::cerr << "Divide by zero warning.";
    // std::exit(0);
  }
  q_next = (1.0/H_OUT(j+1,j)) * q_next;

  /// store q_{j+1} in (j+1)th column of Q
  Q_OUT.setCol(j+1,q_next);

}

void applyGMRES(const Matrix& A, const Vector& b, Vector& x_OUT) {
  // This function executes the GMRES method to solve the square linear system
  // Ax=b.

  // first, we get the size of the system, which should be square
  assert( A.size(0) == A.size(1) );
  const index m{A.size(0)};

  // and we begin the iterative process with
  index n{1};
  // which will specify the dimension of the Krylov subspace we minimize over

  // We declare and initialize the residual vector r = b-Ax (but we assume x=0),
  // our scalar beta = ||r||_2, and the matries Q_{n+1} and H_n and H_n_LSQ
  const double beta{ b.getNorm() };
  Vector r{m};
  r  =   b;
  Matrix Q_np1  { m,  n+1};
  Matrix H_n    {n+1,  n };
  Matrix H_n_LSQ{n+1,  n };

  // for debugging
  std::cout << "\nInitial residual = " << r.getNorm() << '\n';

  // store the normalized residual in first column of Q_{n+1}
  // (a hassle b/c C is row major)
  // this is the first vector in our orthonormal basis for the Krylov subspace
  for (index i{}; i < m; ++i) {
    Q_np1(i,0) = (1.0/r.getNorm())*r(i);
  }

  // now, we iterate over n = 1, 2, 3, ..., n_max until the
  // desired tolerance is acheived
  while (r.getNorm() > consts::epsilon) {

    // for debugging
    std::cout << "\nn = " << n << '\n';

    // apply nth step of Arnoldi Algorithm
    applyArnoldi(A, Q_np1, H_n, n-1);
    H_n_LSQ = H_n;
    // std::cout << "\nQ_np1" << Q_np1 << "\nH_n" << H_n << "\nH_n_LSQ" << H_n_LSQ;

    // find y that minimizes ||beta*e1 - H_n*y||
    Vector v_diag{n};
    Vector beta_e1{n+1};
    beta_e1(0) = beta;
    decompQR(H_n_LSQ, beta_e1, v_diag); // this routine alters H_n!
    // std::cout << "\nH_n" << H_n << "\nH_n_LSQ" << H_n_LSQ;

    Vector y{n};
    solveTriLSQ(H_n_LSQ, beta_e1, y);
    // std::cout << "y" << y;

    // take x_n = x_0 + Q_n*y = Q_n*y
    // we perform the matrix multiplication x = Q_n*y manually since
    // Q_np1 has n+1 many columns whereas y has length n
    for (index i{}; i < m; ++i) {
      x_OUT(i) = 0.0;
      for (index j{}; j < n; ++j) {
        x_OUT(i) += Q_np1(i,j)*y(j);
      }
    }

    // for debugging
    std::cout << "\nx_" << n << " =" << x_OUT;

    // recompute residual after checking for monotone convergence
    if ( (A*x_OUT - b).getNorm() >= r.getNorm() ) {
      std::cout << "\nResidual has increased from " << r.getNorm() << " to "
                << (A*x_OUT - b).getNorm() << ", but the method should "
                << "converge monotonically. Exiting the routine.\n\n";
      std::exit(0);
    }
    else {
      r = A*x_OUT - b;
      std::cout << "\nResidual = " << r.getNorm() << '\n';
      // beta = (A*x_OUT - b).getNorm();
      // std::cout << "\nResidual = " << beta << '\n';
    }

    // check whether we have not reached the theoretical max
    if (n < m) {
      n += 1;
    }
    else {
      std::cerr << "\nMaximum number of iterations reached w/o convergence.";
      std::cerr << "\nExiting the iterative process.\n";
      break;
    }

    // resize Q_{n+1} and H_n for next iteration
    Q_np1.resize(    m , n+1 );
    H_n.resize(     n+1,  n  );
    H_n_LSQ.resize( n+1,  n  );
  }
}
