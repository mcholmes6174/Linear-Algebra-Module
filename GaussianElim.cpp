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
#include <cmath>

void pivotElim(double mat_A[][consts::n], double vec_b[],
               const int m, const int n) {
  // loop over all columns (func. terminates at end of this loop)
  for (int j{}; j<n; ++j) {
    // search for better pivot (p) than current pivot (j)
    int p{j}; // index of pivot row
    for (int k{j+1}; k<m; ++k) {
      if (std::abs(mat_A[k][j]) > std::abs(mat_A[p][j])) p = k;
    }
    // if better pivot (p) is found, then swap rows j and p IN ENTIRE SYSTEM
    if (p != j) tools::swapRows(mat_A,vec_b,j,p);
    // if pivot is still zero, then matrix must be singular
    tools::catchSingular(mat_A[j][j]);
    // if pivot is nonzero, then we loop over rows below row j
    for (int i{j+1}; i<m; ++i) {
      // for each row, we perform the Guassian Elimination row operation
      double scalar{ mat_A[i][j]/mat_A[j][j] };
      for (int k{}; k<n; ++k) {
        mat_A[i][k] -= mat_A[j][k] * scalar;
      }
      // repreat the row operation to vector b to ensure an equivalent system
      vec_b[i] -= vec_b[j] * scalar;
    }
  }
}

void backSub(const double mat_U[][consts::n],
             const double vec_y[], double vec_x[],
             const int m, const int n) {
  // loop over the rows, bottom to top, computing the entries x[i]
  for (int i{m-1}; i>=0; --i) {
    // test for singularity
    tools::catchSingular(mat_U[i][i]);
    // Take the sum of products in row i of Ux=y for j > i to solve for x[i]
    double sum{0.0};
    for (int j{i+1}; j<n; ++j) {
      sum += mat_U[i][j] * vec_x[j];
    }
    vec_x[i] = (vec_y[i] - sum)/mat_U[i][i];
  }
}
