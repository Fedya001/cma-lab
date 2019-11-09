#include "validation.h"
#include "report.h"

int main() {
  if (validation::TestAll<double>()) {
    std::cerr << "All tests passed successfully\n";
  } else {
    std::cerr << "Some tests failed\n";
  }

  report::LatexSolutions("latex/solutions.tex");

  return 0;
}