#include "validation.h"

#include <iostream>

int main() {
  if (validation::TestAll<double>()) {
    std::cerr << "All tests passed successfully\n";
  } else {
    std::cerr << "Some tests failed\n";
  }

  return 0;
}