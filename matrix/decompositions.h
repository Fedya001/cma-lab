#pragma once

#include "square_matrix.h"

#include <cstdint>
#include <vector>

template<class T>
struct LDLTDecomposition {
  std::vector<int8_t> diagonal;
  SquareMatrix<T> low;

  void MultiplyDiagonal();
};

template<class T>
void LDLTDecomposition<T>::MultiplyDiagonal() {
  size_t dim = low.GetDim();
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      low[row][column] *= diagonal[column];
    }
  }
}

template<class T>
struct DLUDecomposition {
  std::vector<size_t> rows_permutations;
  SquareMatrix<T> low_up;
};