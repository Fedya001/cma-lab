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
  int32_t dim = low.GetDim();
  for (int32_t row = 0; row < dim; ++row) {
    for (int32_t column = 0; column < dim; ++column) {
      low[row][column] *= diagonal[column];
    }
  }
}

template<class T>
struct DLUDecomposition {
  std::vector<int32_t> rows_permutations;
  SquareMatrix<T> low_up;
};