#include "square_matrix_manager.h"

#include <algorithm>
#include <numeric>

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
typename SquareMatrixManager<T>::LDLTDecomposition SquareMatrixManager<T>::PerformLDLT() const {
  std::vector<bool> diagonal;
  SquareMatrix<T> matrix(5);

  return {diagonal, matrix};
}

template<class T>
typename SquareMatrixManager<T>::DLUDecomposition SquareMatrixManager<T>::PerformDLU(bool swap_rows) const {
  SquareMatrix<T> low_up = matrix_;
  std::vector<std::pair<size_t, size_t>> swaps;

  // TODO: measure time from this point (@flipper)
  size_t dim = matrix_.GetDim();
  for (size_t iteration = 0; iteration < dim; ++iteration) {
    size_t max_index = iteration;

    if (swap_rows) {
      // Search for max element
      for (size_t index = iteration + 1; index < dim; ++index) {
        if (std::abs(low_up.at(max_index).at(iteration)) < std::abs(low_up.at(index).at(iteration))) {
          max_index = index;
        }
      }

      if (max_index != iteration) {
        // perform swap
        swaps.emplace_back(max_index, iteration);
        low_up.SwapRows(iteration, max_index);
      }
    }

    // base element must be non-zero non-zero
    if (matrix_.at(iteration).at(iteration) == 0) {
      std::string message = "Unable to perform DLU decomposition.";
      if (!swap_rows) {
        message += " Consider set swap_rows to true.";
      }
      throw std::logic_error(message);
    }

    // perform gauss step
    // 1. Low matrix
    for (size_t row_index = iteration + 1; row_index < dim; ++row_index) {
      low_up[row_index][iteration] /= low_up[iteration][iteration];
    }

    // 2. Up matrix
    for (size_t i = iteration + 1; i < dim; ++i) {
      for (size_t j = iteration + 1; j < dim; ++j) {
        low_up[i][j] -= low_up[iteration][j] * low_up[i][iteration];
      }
    }
  }

  std::vector<size_t> rows_permutations;
  std::iota(rows_permutations.begin(), rows_permutations.end(), 1);

  std::reverse(rows_permutations.begin(), rows_permutations.end());
  for (const auto& [lhs, rhs] : swaps) {
    std::swap(rows_permutations[lhs], rows_permutations[rhs]);
  }

  return {rows_permutations, low_up};
}

template<class T>
std::vector<T> SquareMatrixManager<T>::SolveSystem(const std::vector<T>& result) const {
  return std::vector<T>(result);
}

// explicit instantiation of templates
template class SquareMatrixManager<double>;
