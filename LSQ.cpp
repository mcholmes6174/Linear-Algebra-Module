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

 using std::size_t;

 void decompQR(mat_t& A_OUT, mat_t& Q_OUT) {
   // This function computes the full QR decomposition (A=QR) of an m by n
   // (m > n) matrix A using the method of Householder transformations. in the
   // decomposition, Q is m by m and orthogonal, while R is m by n and upper-
   // triangular.
   //
   // The matrix A is overwritten with the matrix R, while the matrix Q is filled
   // with the orthonormalized vectors q_1 through q_m on exit.
 }

 void solveTriLSQ(const mat_t& Q, const mat_t& R, const vec_t& b, vec_t& x_OUT) {
   // This function solves the triangular least squares problem (QR)x = b, where
   // QR is the full QR decomposition of an m by n (m > n) matrix A. It performs
   // backsubstitution to solve the triangular square system R_t * x = Q^T * b,
   // where R_t is the first n rows of R and Q^T is the transpose of Q. It writes
   // the solution into the vector x_OUT on exit.
 }
