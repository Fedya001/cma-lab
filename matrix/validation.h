#pragma once

#include "decompositions.h"
#include "load.h"
#include "matrix_factory.h"
#include "square_matrix.h"

#include <cmath>
#include <iostream>

namespace validation {

// C++17
inline const double EPSILON = 1e-6;

template<class Function>
bool CheckLogicErrorExceptionThrown(Function function) {
  try {
    function();
  } catch (const std::logic_error& ex) {
    // everything is correct
    return true;
  }
  return false;
}

template<class T>
void ValidateDlU(const DLUDecomposition<T>& decomposition,
                 const SquareMatrix<T>& matrix) {
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

      if (std::abs(matrix[inverse_permutation[row]][column] - sum) > EPSILON) {
        throw std::runtime_error("Invalid DLU decomposition\n");
      }
    }
  }
}

template<class T>
void ValidateLDLT(const LDLTDecomposition<T>& decomposition,
                  const SquareMatrix<T>& matrix) {
  size_t dim = matrix.GetDim();
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      T sum = T();
      for (size_t index = 0; index <= std::min(row, column); ++index) {
        sum += decomposition.low[row][index] * decomposition.low[column][index]
            * decomposition.diagonal[index];
      }

      if (std::abs(matrix[row][column] - sum) > EPSILON) {
        throw std::runtime_error("Invalid LDLT decomposition\n");
      }
    }
  }
}

template<class T>
void ValidateLowdiagInverse(const SquareMatrix<T>& inverse,
                  const SquareMatrix<T>& matrix) {
  const auto result = inverse * matrix; // identity expected

  size_t dim = matrix.GetDim();
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      T diff = result[row][column];
      if (row == column) {
        diff -= T(1);
      }
      if (std::abs(diff) > EPSILON) {
        throw std::runtime_error("Invalid lowdiag inverse matrix\n");
      }
    }
  }
}

template<class T>
bool TestAll() {
  // Step out of a cmake-build-debug directory
  const std::string dlu_directory("../data/tests/DLU");
  const std::string ldlt_directory("../data/tests/LDLT");

  const std::string lowdiag_task_directory("../data/task1");

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

    // 4. Test Lowdiag
    {
      auto non_degenerate = loader::LoadMatrix<T>(lowdiag_task_directory + "/sampleA.data");
      manager.SetMatrix(non_degenerate);
      ValidateLowdiagInverse(manager.InverseLowdiag(true), non_degenerate);

      auto degenerate = loader::LoadMatrix<T>(lowdiag_task_directory + "/sampleB.data");
      manager.SetMatrix(degenerate);
      if (!CheckLogicErrorExceptionThrown([&]() {manager.InverseLowdiag(true);})) {
        throw std::runtime_error("Inverse lowdiag fail: can't inverse degenerate matrix\n");
      }
    }

    for (size_t dim : {10, 20, 50}) {
      auto matrix = matrix_factory.CreateLowdiagMatrix(dim);
      manager.SetMatrix(matrix);
      ValidateLowdiagInverse(manager.InverseLowdiag(true), matrix);
    }
  } catch (std::exception& ex) {
    std::cerr << ex.what() << std::endl;
    return false;
  }

  return true;
}

} // namespace validation

