#include <string>

#ifndef IN_OUT
#define IN_OUT

namespace inout {

  using std::size_t;

  std::string askFileChoice(const char);

  int askMethodChoice(const size_t, const size_t);

  void generateSymmSystem(const size_t, const size_t);

}

#endif
