/* Written by Matthew Holmes
 * December 2021
 *
 * This program contains the two necessary routines for performing Guassian
 * Elimination to solve a linear system Ax=b:
 *
 * 1) Gaussian Elimination with Partial Pivoting in the original system Ax=b
 * 2) Back-Substitution in the transformed system Ux=y
 *
 */
#include "GaussianElim.h"
#include "LinAlgToolkit.h"
#include "Matrix.h"
#include "Types.h"
#include "Vector.h"
#include <cmath>

void pivotElim(Matrix& A_OUT, Vector& b_OUT) {

  // get size of system
  index m{ A_OUT.size(0) };
  index n{ A_OUT.size(1) };

  // loop over all columns (func. terminates at end of this loop)
  for (index j{}; j < n; ++j) {

    // search for better pivot (p) than current pivot (j)
    index p{j}; // index of pivot row
    for (auto k{j+1}; k < m; ++k) {
      if (std::abs(A_OUT(k,j)) > std::abs(A_OUT(p,j))) p = k;
    }

    // if better pivot (p) is found, then swap rows j and p IN ENTIRE SYSTEM
    if (p != j) tlk::swapRows(A_OUT,b_OUT,j,p);

    // if pivot is still zero, then matrix must be singular
    tlk::catchSingular(A_OUT(j,j));

    // if pivot is nonzero, then we loop over rows below row j
    for (auto i{j+1}; i<m; ++i) {

      // for each row, we perform the Guassian Elimination row operation
      double scalar{ A_OUT(i,j)/A_OUT(j,j) };
      for (index k{}; k < n; ++k) {
        A_OUT(i,k) -= A_OUT(j,k)*scalar;
      }

      // repreat the row operation to vector b to ensure an equivalent system
      b_OUT(i) -= b_OUT(j)*scalar;

    }
  }
}

void backSub(const Matrix& U, const Vector& y, Vector& x_OUT) {

  // get size of system
  index m{ U.size(0) };
  index n{ U.size(1) };

  // loop over the rows, bottom to top, computing the entries x[i]
  for (auto i{m}; i-- > 0; ) { // loop takes odd shape due to unsigned integers

    // test for singularity
    tlk::catchSingular(U(i,i));

    // Take the sum of products in row i of Ux=y for j > i to solve for x[i]
    double sum{0.0};
    for (auto j{i+1}; j < n; ++j) {
      sum += U(i,j) * x_OUT(j);
    }
    x_OUT(i) = (y(i)-sum)/U(i,i);

  }
}
