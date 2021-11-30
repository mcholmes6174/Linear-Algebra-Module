// namespace to hold important constant values
#ifndef CONSTANTS
#define CONSTANTS

namespace consts{
  inline constexpr double epsilon{5*2.22e-16}; // IEEE machine precision *5
  inline constexpr double tolerance{1.0e-8};   // stopping condition
  inline constexpr int m{4}; // num rows
  inline constexpr int n{4}; // num cols
}

#endif
