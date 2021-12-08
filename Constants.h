// namespace to hold important constant values
#ifndef CONSTANTS
#define CONSTANTS

namespace consts{
  inline constexpr double epsilon{5*2.22e-16}; // IEEE machine precision *5
  inline constexpr double tolerance{1.0e-8};   // stopping condition
  inline constexpr int    dispPrec{8};         // precision to display values
  inline constexpr int    writePrec{16};       // precision to write values
}

#endif
