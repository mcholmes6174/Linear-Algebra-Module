#include <string>

#ifndef IN_OUT
#define IN_OUT

namespace inout {

  std::string askFileChoice(const char);

  int askMethodChoice();

  void generateSymmSystem(const size_t, const size_t);

}

#endif
