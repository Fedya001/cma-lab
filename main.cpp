#include "square_matrix.h"
#include "square_matrix_manager.h"
#include "print_utils.h"
#include "validation.h"

#include <iostream>

int main() {
  SquareMatrix<double> matrixA({
                                  {4, 2, 1},
                                  {5, 10, 3},
                                  {4, 5, 4}
                              });

  auto decompositionA = SquareMatrixManager(matrixA).PerformDLU(true);

  PrintMatrix(std::cout, decompositionA.low_up, false);
  std::cout << std::endl;
  PrintColumn(std::cout, decompositionA.rows_permutations, false);

  SquareMatrix<double> matrixB({
                                            {4, 5, 7},
                                            {5, 1, 2},
                                            {7, 2, 10}
                                        });
  auto decompositionB = SquareMatrixManager(matrixB).PerformLDLT();
  PrintMatrix(std::cout, decompositionB.low, false, 4);
  std::cout << std::endl;
  PrintColumn(std::cout, decompositionB.diagonal);
  std::cout << std::endl;

  try {
    validation::ValidateDlU(decompositionA, matrixA);
    validation::ValidateLDLT(decompositionB, matrixB);
  } catch (std::exception& error) {
    std::cerr << error.what() << std::endl;
  }

  if (validation::TestAll<double>()) {
    std::cerr << "All tests passed successfully\n";
  } else {
    std::cerr << "Some tests failed\n";
  }

  // test permutations
  matrixA.SwapRows({1, 2, 0});
  PrintMatrix(std::cout, matrixA);

  SquareMatrix<double> test({
                                   {1, 1, 1, 1, 1},
                                   {2, 2, 2, 2, 2},
                                   {3, 3, 3, 3, 3},
                                   {4, 4, 4, 4, 4},
                                   {5, 5, 5, 5, 5}
                               });
  test.SwapRows({1, 0, 3, 4, 2});
  PrintMatrix(std::cout, test);

  return 0;
}