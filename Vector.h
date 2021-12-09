// This file contains the declaration of the Matrix class.
#include "Types.h"

#ifndef VECTOR
#define VECTOR

// here we define a new class called Vector
class Vector {

// private members can only be accessed by other members of Vector
private:
  // we create two private member variables of types std::size_t and std::vector
  index m_row{};
  vec_t m_vec{};

// public members are accessible by anyone; we create many public member funcs
public:

  // here is a default constructor able to accept one user-provided value
  // **we use member initializer lists to initialize our member variables
  //   instead of assignment within the function body**
  Vector(const index m=1) : m_row{m}, m_vec(m) {
    // when a new Vector object is created, this default constructor is called
    // and initializes m_vec with length m.
  }

  // trivial member function to return the size of the vector
  index  getSize() { return m_vec.size(); }

  // the following member functions are defined within the Vector.cpp file
  // because they are > 1 line of code

  void   setSize(const index m);

  void   setVal(const index i, const double v_i);

  double getVal(const index i);

  void   show();

  double getNorm(const double p = 2.0);
};

#endif
