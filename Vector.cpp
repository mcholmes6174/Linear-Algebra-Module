// This file contains the implementations of the member functions for the
// Vector class
#include "Constants.h"
#include "LinAlgToolkit.h"
#include "Matrix.h"
#include "Types.h"
#include "Vector.h"
#include <cassert>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>

void Vector::resize(const index m) {
  // to reset the length of the vector
  m_row = m;
  m_vec.resize(m_row);
}

void Vector::normalize() {
  double norm{ getNorm() };
  for (index i{}; i < m_vec.size(); ++i) {
    m_vec[i] /= norm;
  }
}

void Vector::show() const {
  // to print the vector to the screen
  using std::cout;
  cout << std::setprecision(consts::dispPrec);
  cout << '\n';
  for (auto v_i : m_vec) {
    cout << '\n' << std::setw(consts::dispPrec+2)
         << std::right << std::showpoint << v_i;
  }
  cout << '\n';
}

double Vector::getNorm(const double p) const {
  // to return the p-norm (2 by default) norm of the vector
  double norm{};
  for (auto v_i : m_vec) {
    norm += std::pow(v_i, p);
  }
  return std::pow(norm, 1.0/p);
}

void Vector::load(const std::string filename) {
  // this function reads the contents of "filename" and stores the data in
  // the vector x_OUT.
  // we'll read from a file called "filename"
  std::ifstream inf{filename};
  // in case we cannot open the input file stream
  tlk::catchFileError(filename, !inf);
  // get length of vector from file
  index m{};
  std::string strVal;
  inf >> strVal;
  m = static_cast<index>(std::stoi(strVal));
  // check that size of x_OUT matches size of vector in file
  if (m != m_vec.size()) {
    std::cout << "\nError: Size of vector x_OUT in readVector() does not "
              << "match size of vector in " << filename << ".\n";
    std::exit(0);
  }
  // loop through file inserting values into vector
  for (index i{}; i < m; ++i) {
    // get value and store within A
    inf >> strVal;
    m_vec[i] = std::stod(strVal);
  }
}

void Vector::write(const std::string filename) const {
  // this function writes the contents of a vector x into "filename"
  // create a file
  std::ofstream outf{filename};
  // in case we cannot open the output file stream
  if (!outf) {
    std::cerr << "\nError: soln.txt could not be opened for writing\n\n";
    return;
  }
  // get length of vector and write into file
  index m{m_vec.size()};
  outf << m << '\n';
  // write into the file
  std::setprecision(consts::writePrec);
  for (auto v_i : m_vec) {
    outf << std::setw(consts::writePrec+2)
         << std::left << std::showpoint << v_i << '\n';
  }
  // when outf goes out of scope, the ofstream destructor will close the file
}

//*****************************************************************************/
// here we define member functions in order to perform operator overloading
//*****************************************************************************/

double& Vector::operator()(const index i) {
  // here we allow subscripting on our vector objects
  assert(i < m_row);
  return m_vec[i];
}

double  Vector::operator()(const index i) const {
  // here we allow subscripting on our *constant* vector objects
  assert(i < m_row);
  return m_vec[i];
}

Vector& Vector::operator=(const Vector& x) {
  // here we create an assignment operator
  if (this == &x) return *this; // self-assignment guard

  m_row = x.m_row;
  m_vec = x.m_vec;

  return *this;
}

//*****************************************************************************/
// below we define *normal* functions in order to perform operator overloading
//*****************************************************************************/

Vector operator-(const Vector& x) {
  // here we negate all the elements of a vector x
  Vector minus_x{x.size()};
  for (index i{}; i < x.size(); ++i) {
    minus_x(i) = -x(i);
  }
  return minus_x;
}

Vector operator+(const Vector& x, const Vector& y) {
  // here we take the componentwise sum of two vectors x and y
  assert(x.size() == y.size());
  Vector z{x.size()};
  for (index i{}; i < x.size(); ++i) {
    z(i) = x(i) + y(i);
  }
  return z;
}

Vector operator-(const Vector& x, const Vector& y) {
  // here we take the componentwise difference of two vectors x and y
  assert(x.size() == y.size());
  return x + (-y); // we can now use the Vector operators defined above
}

Vector operator*(const double scalar, const Vector& x) {
  // here we scale a vector
  Vector scaled_x{x.size()};
  for (index i{}; i < x.size(); ++i) {
    scaled_x(i) = scalar*x(i);
  }
  return scaled_x;
}

Vector operator*(Vector& x, const double scalar) {
  // we want the scaling operation to be symmetric
  return scalar*x;
}
