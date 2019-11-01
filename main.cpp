#include "square_matrix.h"
#include "square_matrix_manager.h"
#include "print_utils.h"
#include "validator.h"

#include <iostream>

int main() {
  SquareMatrix<double> matrix({
                                  {4, 2, 1},
                                  {5, 10, 3},
                                  {4, 5, 4}
                              });

  auto result = SquareMatrixManager(matrix).PerformDLU(true);
  PrintMatrix(std::cout, result.low_up, false);
  std::cout << std::endl;
  PrintColumn(std::cout, result.rows_permutations, false);

  SquareMatrix<double> symmetric_matrix({
                                            {4, 5, 7},
                                            {5, 1, 2},
                                            {7, 2, 10}
                                        });
  auto decomposition = SquareMatrixManager(symmetric_matrix).PerformLDLT();
  PrintMatrix(std::cout, decomposition.low, false, 4);
  std::cout << std::endl;
  PrintColumn(std::cout, decomposition.diagonal);

  return 0;
}