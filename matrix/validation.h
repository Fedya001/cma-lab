#pragma once

#include "decompositions.h"
#include "load.h"
#include "matrix_factory.h"
#include "square_matrix.h"

#include <cmath>
#include <iostream>

namespace validation {

template<class T>
void ValidateDlU(const DLUDecomposition<T>& decomposition,
                 const SquareMatrix<T>& matrix) {
  const double epsilon = 1e-6;
  size_t dim = matrix.GetDim();

  std::vector<size_t> inverse_permutation(dim);
  for (size_t index = 0; index < dim; ++index) {
    inverse_permutation[decomposition.rows_permutations[index]] = index;
  }

  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      T sum = T();
      for (size_t index = 0; index <= std::min(row, column); ++index) {
        auto low_element = decomposition.low_up[row][index];
        if (index == row) {
          low_element = 1;
        }
        sum += low_element * decomposition.low_up[index][column];
      }

      if (std::abs(matrix[inverse_permutation[row]][column] - sum) > epsilon) {
        throw std::runtime_error("Invalid DLU decomposition\n");
      }
    }
  }
}

template<class T>
void ValidateLDLT(const LDLTDecomposition<T>& decomposition,
                  const SquareMatrix<T>& matrix) {
  const double epsilon = 1e-6;
  size_t dim = matrix.GetDim();
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      T sum = T();
      for (size_t index = 0; index <= std::min(row, column); ++index) {
        sum += decomposition.low[row][index] * decomposition.low[column][index]
            * decomposition.diagonal[index];
      }

      if (std::abs(matrix[row][column] - sum) > epsilon) {
        throw std::runtime_error("Invalid LDLT decomposition\n");
      }
    }
  }
}

template<class T>
bool TestAll() {
  // Step out of a cmake-build-debug directory
  const std::string dlu_directory("../data/tests/DLU");
  const std::string ldlt_directory("../data/tests/LDLT");

  auto manager = SquareMatrixManager(SquareMatrix<T>(0));

  try {
    // 1. Test DLU
    for (const auto& matrix : loader::LoadMatrices<T>(dlu_directory)) {
      manager.SetMatrix(matrix);
      ValidateDlU(manager.PerformDLU(true), matrix);
    }

    // 2. Test LDLT
    for (const auto& matrix : loader::LoadMatrices<T>(ldlt_directory)) {
      manager.SetMatrix(matrix);
      ValidateLDLT(manager.PerformLDLT(), matrix);
    }

    // 3. // Random tests on DLU and LDLT
    MatrixFactory<T> matrix_factory(-1000, 1000);
    for (size_t dim : {10, 20, 50, 100}) {
      auto matrix = matrix_factory.CreateRandomMatrix(dim, false);
      manager.SetMatrix(matrix);
      ValidateDlU(manager.PerformDLU(true), matrix);

      auto symmetric_matrix = matrix_factory.CreateRandomMatrix(dim, true);
      manager.SetMatrix(symmetric_matrix);
      ValidateLDLT(manager.PerformLDLT(), symmetric_matrix);
    }
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return false;
  }

  return true;
}

} // namespace validation

