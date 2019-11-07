#include "relaxation.h"
#include "validation.h"
#include "print_utils.h"

#include <iostream>

int main() {
  if (validation::TestAll<double>()) {
    std::cerr << "All tests passed successfully\n";
  } else {
    std::cerr << "Some tests failed\n";
  }

  PrintRow(std::cout, SolveSystem(5, 0.1));

  return 0;
}