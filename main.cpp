#include "square_matrix.h"
#include "square_matrix_manager.h"

#include <iostream>

int main() {
  SquareMatrix<double> matrix({
                                  {4, 2, 1},
                                  {5, 10, 3},
                                  {4, 5, 4}
                              });

  auto result = SquareMatrixManager(matrix).PerformDLU();
  PrintMatrix(std::cout, result.low_up, true);

  return 0;
}