// This file contains the implementations of the member functions for the
// Vector class
#include "Constants.h"
#include "Types.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
#include <iomanip>
#include <iostream>

void Vector::setSize(const index m) {
  // to reset the length of the vector
  m_row = m;
  m_vec.resize(m_row);
}

void Vector::setVal(const index i, const double v_i) {
  // to set the value of a particular entry
  assert(i < m_row);
  m_vec[i] = v_i;
}

double Vector::getVal(const index i) {
  // to get the value stored at index i
  assert(i < m_row);
  return m_vec[i];
}

void Vector::show() {
  // to print the vector to the screen
  using std::cout;
  cout << std::setprecision(consts::dispPrec);
  for (auto v_i : m_vec) {
    cout << '\n' << std::setw(consts::dispPrec+2)
         << std::right << std::showpoint << v_i;
  }
  cout << '\n';
}

double Vector::getNorm(const double p) {
  // to return the p-norm (2 by default) norm of the vector
  double norm{};
  for (auto v_i : m_vec) {
    norm += std::pow(v_i, p);
  }
  return    std::pow(norm, 1.0/p);
}
