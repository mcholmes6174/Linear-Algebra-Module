/* Written by Matthew Holmes
 * November 2021
 *
 * This program contains the two necessary routines for performing Guassian
 * Elimination to solve a linear system Ax=b:
 *
 * 1) Gaussian Elimination with Partial Pivoting in the original system Ax=b
 * 2) Back-Substitution in the transformed system Ux=y
 *
 */
#include "GaussianElim.h"
#include "LinAlgTools.h"
#include "Types.h"
#include <cmath>
#include <vector>

void pivotElim(mat_t& A, vec_t& b) {
  // get size of system
  int m{std::size(A)};
  int n{std::size(A[0])};
  // loop over all columns (func. terminates at end of this loop)
  for (int j{}; j<n; ++j) {
    // search for better pivot (p) than current pivot (j)
    int p{j}; // index of pivot row
    for (int k{j+1}; k<m; ++k) {
      if (std::abs(A[k][j]) > std::abs(A[p][j])) p = k;
    }
    // if better pivot (p) is found, then swap rows j and p IN ENTIRE SYSTEM
    if (p != j) tools::swapRows(A,b,j,p);
    // if pivot is still zero, then matrix must be singular
    tools::catchSingular(A[j][j]);
    // if pivot is nonzero, then we loop over rows below row j
    for (int i{j+1}; i<m; ++i) {
      // for each row, we perform the Guassian Elimination row operation
      double scalar{ A[i][j]/A[j][j] };
      for (int k{}; k<n; ++k) {
        A[i][k] -= A[j][k] * scalar;
      }
      // repreat the row operation to vector b to ensure an equivalent system
      b[i] -= b[j] * scalar;
    }
  }
}

void backSub(const mat_t& U, const vec_t& y, vec_t& x) {
  // get size of system
  int m{std::size(U)};
  int n{std::size(U[0])};
  // loop over the rows, bottom to top, computing the entries x[i]
  for (int i{m-1}; i>=0; --i) {
    // test for singularity
    tools::catchSingular(U[i][i]);
    // Take the sum of products in row i of Ux=y for j > i to solve for x[i]
    double sum{0.0};
    for (int j{i+1}; j<n; ++j) {
      sum += U[i][j] * x[j];
    }
    x[i] = (y[i] - sum)/U[i][i];
  }
}
