#pragma once

#include "square_matrix.h"
#include "decompositions.h"
#include "load_matrix.h"
#include "matrix_generator.h"

#include <cmath>
#include <experimental/filesystem>
#include <iostream>

namespace validation {

namespace filesystem = std::experimental::filesystem;

template<class T>
void ValidateDlU(const DLUDecomposition<T>& decomposition,
                 const SquareMatrix<T>& matrix) {
  const double epsilon = 1e-6;
  size_t dim = matrix.GetDim();

  std::vector<size_t> inverse_permutation(dim);
  for (size_t index = 0; index < dim; ++index) {
    inverse_permutation[decomposition.rows_permutations.at(index)] = index;
  }

  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      T sum = T();
      for (size_t index = 0; index <= std::min(row, column); ++index) {
        auto low_element = decomposition.low_up.at(row).at(index);
        if (index == row) {
          low_element = 1;
        }
        sum += low_element * decomposition.low_up.at(index).at(column);
      }

      if (std::abs(matrix.at(inverse_permutation.at(row)).at(column) - sum) > epsilon) {
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
      for (size_t index = 0; index < std::min(row, column); ++index) {
        sum += decomposition.low.at(row).at(index) * decomposition.low.at(column).at(index);
      }

      if (std::abs(matrix.at(row).at(column) - sum) > epsilon) {
        throw std::runtime_error("Invalid LDLT decomposition\n");
      }
    }
  }
}

template<class T>
bool TestMatrix(const SquareMatrix<T>& matrix) {
  auto manager = SquareMatrixManager(matrix);

  auto dlu = manager.PerformDLU(true);
  auto ldlt = manager.PerformLDLT();

  try {
    ValidateDlU(dlu, matrix);
    ValidateLDLT(ldlt, matrix);
  } catch (std::runtime_error& ex) {
    std::cerr << ex.what();
    return false;
  }

  return true;
}

template<class T>
bool TestAll() {
  const std::string directory("data/tests");

  std::vector<std::string> test_files;
  try {
    if (filesystem::exists(directory) && filesystem::is_directory(directory)) {
      filesystem::recursive_directory_iterator iter(directory);

      // Iterate till end
      while (iter != filesystem::recursive_directory_iterator()) {
        if (filesystem::is_directory(iter->path())) {
          // c++17 Filesystem API to skip current directory iteration
          iter.disable_recursion_pending();
        } else {
          // Add the name in vector
          test_files.push_back(iter->path().string());
        }

        std::error_code error;
        iter.increment(error);
        if (error) {
          std::cerr << "Error While Accessing : " << iter->path().string()
                    << " :: " << error.message() << '\n';
        }
      }
    }
  }
  catch (std::system_error& e) {
    std::cerr << e.what();
  }

  // Tests in data/tests
  for (const auto& test_file : test_files) {
    if (!TestMatrix(loader::LoadMatrix<T>(test_file))) {
      return false;
    }
  }

  // Random tests
  MatrixGenerator<T> matrix_generator(-1000, 1000);
  for (size_t dim : {10, 20, 50, 200, 500}) {
    if (!TestMatrix(matrix_generator.generate(dim))) {
      return false;
    }
  }

  return true;
}

} // namespace validation

