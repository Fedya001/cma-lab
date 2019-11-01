#pragma once

#include "square_matrix.h"

#include <fstream>

namespace loader {

template<class T>
SquareMatrix<T> LoadMatrix(const std::string& filename) {
  std::ifstream input(filename);

  size_t dim;
  input >> dim;

  T element;
  SquareMatrix<T> matrix(dim);
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      input >> element;
      matrix[row][column] = element;
    }
  }

  return matrix;
}

} // namespace loader
