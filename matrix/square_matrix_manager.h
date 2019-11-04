#pragma once

#include "decompositions.h"
#include "square_matrix.h"
#include "matrix_factory.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <functional>

template<class T>
class SquareMatrixManager {
 public:
  enum class MatrixType {
    UPPER_TRIANGULAR,
    LOWER_TRIANGULAR,
    LOWDIAG,
    SYMMETRIC,
    ORDINARY
  };

  explicit SquareMatrixManager(SquareMatrix<T> matrix);

  SquareMatrix<T> GetMatrix() const;
  void SetMatrix(const SquareMatrix<T>& matrix);

  LDLTDecomposition<T> PerformLDLT() const;
  DLUDecomposition<T> PerformDLU(bool swap_rows = false) const;

  // Lowdiag = Lower triangular matrix + one extra diagonal above it
  SquareMatrix<T> InverseLowdiag(bool swap_rows = false) const;

  std::vector<T> SolveSystem(const std::vector<T>& result) const;

 private:
  SquareMatrix<T> matrix_;

  void EnsureMatrixType(MatrixType matrix_type) const;
 };

template<class T>
SquareMatrixManager<T>::SquareMatrixManager(SquareMatrix<T> matrix)
    : matrix_(std::move(matrix)) {}

template<class T>
SquareMatrix<T> SquareMatrixManager<T>::GetMatrix() const {
  return matrix_;
}

template<class T>
void SquareMatrixManager<T>::SetMatrix(const SquareMatrix<T>& matrix) {
  matrix_ = matrix;
}

template<class T>
LDLTDecomposition<T> SquareMatrixManager<T>::PerformLDLT() const {
  EnsureMatrixType(MatrixType::SYMMETRIC);

  size_t dim = matrix_.GetDim();
  std::vector<int8_t> diagonal(dim, 1);
  SquareMatrix<T> matrix = matrix_;

  // We store the resulting matrix in transposed vector with necessary zeros.
  // This operations take O(n^2) while the whole algorithm
  // takes O(n^3 / 2) that is why there is no problem
  std::vector<std::vector<T>> transposed;

  for (size_t iteration = 0; iteration < dim; ++iteration) {
    if (matrix[iteration][iteration] == 0) {
      throw std::logic_error("Unable to perform LDLT decomposition.");
    }

    for (size_t row = iteration + 1; row < dim; ++row) {
      auto coefficient = -matrix[iteration][row] / matrix[iteration][iteration];
      matrix.AddRowToOther(row, iteration, coefficient, row, dim - 1);
    }

    // Form a row of resulting matrix low
    std::vector<T> row;
    for (size_t index = 0; index <= iteration; ++index) {
      row.push_back(matrix[index][iteration]);
    }
    row.resize(dim, T(0));
    transposed.push_back(row);
  }

  // Finish LDLT: Divide on square root
  for (size_t column = 0; column < dim; ++column) {
    if (transposed[column][column] < 0) {
      diagonal[column] = -1;
    }
    T root = sqrt(diagonal[column] * matrix[column][column]);

    for (size_t row = 0; row < dim; ++row) {
      transposed[row][column] /= diagonal[column] * root;
    }
  }

  return {diagonal, SquareMatrix<T>(std::move(transposed))};
}

template<class T>
DLUDecomposition<T> SquareMatrixManager<T>::PerformDLU(bool swap_rows) const {
  SquareMatrix<T> low_up = matrix_;
  std::vector<std::pair<size_t, size_t>> swaps;

  // TODO: measure time from this point (@flipper)
  size_t dim = matrix_.GetDim();
  for (size_t iteration = 0; iteration < dim; ++iteration) {

    if (swap_rows) {
      size_t max_index = iteration;
      // Search for max element
      for (size_t index = iteration + 1; index < dim; ++index) {
        if (std::abs(low_up[max_index][iteration]) < std::abs(low_up[index][iteration])) {
          max_index = index;
        }
      }

      if (max_index != iteration) {
        swaps.emplace_back(max_index, iteration);
        low_up.SwapRows(max_index, iteration);
      }
    }

    // Base element must be non-zero
    if (low_up[iteration][iteration] == 0) {
      std::string message = "Unable to perform DLU decomposition.";
      if (!swap_rows) {
        message += " Consider set swap_rows to true.";
      }
      throw std::logic_error(message);
    }

    if (iteration != dim - 1) {
      // Perform gauss step
      // 1. Low matrix
      low_up.MultiplyColumn(iteration, T(1) / low_up[iteration][iteration],
                            iteration + 1, dim - 1);


      // 2. Up matrix
      for (size_t row = iteration + 1; row < dim; ++row) {
        low_up.AddRowToOther(row, iteration, -low_up[row][iteration],
                             iteration + 1, dim - 1);
      }
    }
  }

  std::vector<size_t> rows_permutations(dim);
  std::iota(rows_permutations.begin(), rows_permutations.end(), 0);

  std::reverse(swaps.begin(), swaps.end());
  for (const auto&[lhs, rhs] : swaps) {
    std::swap(rows_permutations[lhs], rows_permutations[rhs]);
  }

  return {rows_permutations, low_up};
}

template<class T>
SquareMatrix<T> SquareMatrixManager<T>::InverseLowdiag(bool swap_rows) const {
  EnsureMatrixType(MatrixType::LOWDIAG);
  size_t dim = matrix_.GetDim();

  SquareMatrix<T> left_matrix = matrix_;
  SquareMatrix<T> inverse = MatrixFactory<T>::CreateIdentity(dim);

  const std::string message = "Unable to find inverse of lowdiag matrix.";
  const double epsilon = 1e-6;

  // Step 1. Go up and make zeros on odd diagonal
  for (size_t iter = 0; iter < dim; ++iter) {
    // Smells bad, I know. Suggest a fix, please
    size_t iteration = dim - iter - 1;

    if (swap_rows) {
      if (iteration != 0 && std::abs(left_matrix[iteration][iteration]) <
          std::abs(left_matrix[iteration - 1][iteration])) {
        left_matrix.SwapRows(iteration - 1, iteration);
        inverse.SwapRows(iteration - 1, iteration);
      }
    }

    // Base element must be non-zero
    if (left_matrix[iteration][iteration] == 0) {
      if (!swap_rows) {
        throw std::logic_error(message + " Consider set swap_rows to true.");
      }
      throw std::logic_error(message);
    }

    inverse.MultiplyRow(iteration, T(1) / left_matrix[iteration][iteration]);
    left_matrix.MultiplyRow(iteration, T(1) / left_matrix[iteration][iteration], 0, iteration);

    // If necessary, zero one element above base one
    if (iteration > 0) {
      inverse.AddRowToOther(iteration - 1, iteration, -left_matrix[iteration - 1][iteration]);
      left_matrix.AddRowToOther(iteration - 1, iteration, -left_matrix[iteration - 1][iteration], 0, iteration);
    }
  }

  // Step 2. Go down and make matrix a diagonal
  // There are already a_{ii} = 1 \forall i \in \overline{0, dim - 1}
  for (size_t iteration = 0; iteration < dim; ++iteration) {
    if (std::abs(left_matrix[iteration][iteration] - T(1)) > epsilon) {
      throw std::logic_error(message);
    }
    for (size_t row = iteration + 1; row < dim; ++row) {
      inverse.AddRowToOther(row, iteration, -left_matrix[row][iteration]);
    }
  }

  return inverse;
}

template<class T>
std::vector<T> SquareMatrixManager<T>::SolveSystem(const std::vector<T>& result) const {
  return std::vector<T>(result);
}

template<class T>
void SquareMatrixManager<T>::EnsureMatrixType(MatrixType matrix_type) const {
  std::function<bool(size_t, size_t)> validate_function;
  std::string message;

  switch (matrix_type) {
    case MatrixType::SYMMETRIC:
      validate_function = [this](size_t i, size_t j) {return matrix_[i][j] == matrix_[j][i];};
      message = "Invalid matrix structure (symmetric matrix required)";
      break;
    case MatrixType::UPPER_TRIANGULAR:
      validate_function = [this](size_t i, size_t j) {return j >= i || matrix_[i][j] == 0;};
      message = "Invalid matrix structure (upper triangular matrix required)";
      break;
    case MatrixType::LOWER_TRIANGULAR:
      validate_function = [this](size_t i, size_t j) {return i >= j || matrix_[i][j] == 0;};
      message = "Invalid matrix structure (lower triangular matrix required)";
      break;
    case MatrixType::LOWDIAG:
      validate_function = [this](size_t i, size_t j) {return i + 1 >= j || matrix_[i][j] == 0;};
      message = "Invalid matrix structure (lowdiag matrix required)";
      break;
    default:
      return;
  }

  size_t dim = matrix_.GetDim();
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      if (!validate_function(row, column)) {
        throw std::invalid_argument(message);
      }
    }
  }
}
