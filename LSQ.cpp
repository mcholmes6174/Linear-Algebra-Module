/* Written by Matthew Holmes
* December 2021
*
* This program contains all the routines necessary to perform QR decomposition
* via Householder transformations to solve a least squares problem Ax\approx b.
*
*/
#include "LinAlgToolkit.h"
#include "LSQ.h"
#include "Types.h"
#include <cmath>
#include <iostream>

using std::size_t;

void decompQR(mat_t& A_OUT, vec_t& b_OUT, vec_t& v_diag_OUT) {
  // This function computes the full QR decomposition (A=QR) of an m by n
  // (m > n) matrix A using the method of Householder transformations. in the
  // decomposition, Q is m by m and orthogonal, while R is m by n and upper-
  // triangular.
  //
  // The upper diagonals of the matrix A areoverwritten with the nonzero
  // entires of R, while the columns of A below its diagonal contain the
  // vectors v_j = (0, 0, ..., a_jj + s_j, a_{j+1}j, ..., a_mj)^T used in the
  // creation of the Householder matrices H = I-2vv^T excluding the jth entry
  // a_jj + s_j, which is stored in the separate vector v_diag_OUT. The action
  // of the matrix Q can be fully reconstructed using the vectors v_j. Finally,
  // the householder transformations are also applied to the vector b so that
  // the initial and final systems are equivalent, and we have HA = Hb, or,
  // Q^T*A = Q^T*b

  // we loop over the columns of A
  size_t n{A_OUT[0].size()};
  for (size_t j{}; j < n; ++j) {

  // we get the signed norm of the jth column for i >= j
  double s_j{ std::copysign(1.0, A_OUT[j][j]) };
  double norm{};

  // we also initialize the vector v_j within the loop for efficiency
  size_t m{A_OUT.size()};
  vec_t v_j(m);

  for (size_t i{j}; i < m; ++i) {
    norm  += std::pow( A_OUT[i][j], 2.0 );
    v_j[i] = A_OUT[i][j];
  }
  s_j *= std::sqrt( norm );

  // make sure to add s_j to jth entry of v_j
  v_j[j] += s_j;

  // now we normalize v_j
  tlk::normalizeVec(v_j);

  // next, we apply the matrix multiplication H*A where H = I - 2vv^T
  // we get A = A - 2(vv^T)*A
  for (size_t k{}; k < n; ++k) {

    // get kth column vector of A
    vec_t a_k(m);
    tlk::nthColumn("grab",A_OUT,a_k,k);

    // get inner product < v_j, a_k >
    // and subtract 2*(vv^T)*A from kth column of A
    double vTa{ tlk::innerProd(v_j,a_k) };
    for (size_t i{}; i < m; ++i) {
      A_OUT[i][k] -= 2*v_j[i]*vTa;
    }

  }

  // apply Householder transformation to vector b as well
  double vTb{ tlk::innerProd(v_j,b_OUT) };
  for (size_t i{}; i < m; ++i) {
    b_OUT[i] -= 2*v_j[i]*vTb;
  }

  // store the values of v_j for i > j in A_OUT to save memory
  for (size_t i{j+1}; i < m; ++i) {
    A_OUT[i][j] = v_j[i];
  }

  // finally, we save the jth entry of v_j in v_diag_OUT
  v_diag_OUT[j] = v_j[j];

  }
}

void solveTriLSQ(const mat_t& R, const vec_t& Qb, vec_t& y_OUT) {
  // This function solves the triangular least squares problem (QR)y = b, via
  // backsubstitution of the equation R*y = Q^T*b where Q and R form the reduced
  // QR factorization for A in the equation Ay=b.
  //
  // The inputs should be exactly the outputs of the above decompQR() function
  // with the additional parameter x_OUT to which the solution vector will be
  // written on exit.

  // get size of the square n x n system R*y = Q^T*b
  size_t n{ R[0].size() };

  // loop over the rows of R, bottom to top, computing the entries y[i]
  for (auto i{n}; i-- > 0; ) { // loop takes odd shape due to unsigned integers

    // test for singularity
    tlk::catchSingular(R[i][i]);

    // take the sum of products in row i of Ry=Qb for j > i
    double sum{0.0};
    for (auto j{i+1}; j < n; ++j) {
      sum += R[i][j] * y_OUT[j];
    }

    // and solve for y[i]
    y_OUT[i] = (Qb[i] - sum)/R[i][i];
  }

}
