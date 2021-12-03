/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains all the routines necessary to perform the Generalized
 * Minimal Residual (GMRES) method for solving a square linear system Ax=b.
 *
 * NOTE: The GMRES method can be applied to ANY system where the matrix A is
 *       square and invertible.
 *
 */
 #include "GMRES.h"
 #include "LinAlgToolkit.h"
 #include "LSQ.h"
 #include "Types.h"
 #include <iostream>
 #include <string>

 using std::size_t;

void applyArnoldi(const mat_t& A, mat_t& Q_OUT, mat_t& H_OUT) {
  // This function applies the Arnoldi Algorithm (w/ modified Gram-Schmidt) in
  // order to compute an orthonormal basis for the Krylov subspace K_n(A,b),
  // where the vector b is implicitly contained in the first column of Q as
  // q_1 = r_0/beta where r_0 = b - Ax_0 and beta = ||r_0||_2.
  //
  // The matrix Q is modified directly, with its columns containing the
  // orthonormal basis on exit. The matrix H is modified directly as well,
  // becoming the n+1 by n Hessenberg version of A on exit.

  // get sizes
  size_t m{    Q_OUT.size() };
  size_t n{ Q_OUT[0].size() };

  // loop over columns of Q
  for (size_t j{}; j < n; ++j) {

    // declare and zero-initialize column vectors
    vec_t q(m);                         // q_{j+1}
    vec_t q_j(m);                       // q_j
    tlk::nthColumn("grab",Q_OUT,q_j,j); // get q_j from Q_OUT
    tlk::matVecMul(A,q_j,q);            // q_{j+1} = A*q_j

    // loop under index j and modify q to ensure mutual orthogonality
    for (size_t i{}; i <= j; ++i) {
      H_OUT[i][j] = tlk::innerProd(q,q_j);       // < q_{j+1}, q_i >
      vec_t q_i(m);
      tlk::nthColumn("grab",Q_OUT,q_i,i);
      q = tlk::updateVector(q,-H_OUT[i][j],q_i); // q_{j+1} += -h_{i,j}*q_i
    }

    // normalize the column vector q
    H_OUT[j+1][j] = tlk::getNorm(q);                      // get norm of q_{i+j}
    if (H_OUT[j+1][j] < consts::epsilon) break;            // if zero, stop
    if (j < n) q = tlk::makeVecCopy(q, 1.0/H_OUT[j+1][j] ); // normalize q_{j+1}
    /// store q_{j+1} in (j+1)th column of Q
    tlk::nthColumn("give",Q_OUT,q,j+1);
  }
}

void applyGMRES(mat_t& A_OUT, const vec_t& b, vec_t& x_OUT, const size_t n) {
  // This function executes the GMRES method to solve the square linear system
  // Ax=b. The parameter n specifies the dimension of the Krylov subspace over
  // which to minimize ||Ax-b||_2 in order to obtain a solution x.

  // first, we get the size of the system
  const size_t m{A_OUT.size()};

  // and declare and initialize vector r_0, scalar beta, and matrix U_n
  vec_t  r(m);
  r = tlk::makeVecCopy(b);
  double beta{ tlk::getNorm(r) };
  mat_t  U_n(m,vec_t(n));

  // store (1/beta)*r in first column of U_n (a hassle b/c C is row major)
  for (size_t i{}; auto r_i : r) {
    U_n[i][0] = (1.0/beta)*r_i;
    ++i;
  }

  // for debugging
  std::cout << "r\n";
  tlk::showVector(r);
  std::cout << "U\n";
  tlk::showMatrix(U_n);

  // next, apply the Arnoldi Algorithm to generate
  // an orthogonal basis for the Krylov subspace
  mat_t H_n(n+1,vec_t(n));
  applyArnoldi(A_OUT, U_n, H_n);

  // for debugging
  std::cout << "r\n";
  tlk::showVector(r);
  std::cout << "U\n";
  tlk::showMatrix(U_n);
  std::cout << "H\n";
  tlk::showMatrix(H_n);

  // now we can obtain a QR decomposition for the matrix H_n to solve the
  // least squares problem min_y ||beta*e_1 - H_n*y||_2
  mat_t Q(m,vec_t(n));
  decompQR(H_n,Q);

  // finally, we take the QR decomposition and solve the corresponding
  // triangular least squares problem
  solveTriLSQ(Q,H_n,b,x_OUT);
}
