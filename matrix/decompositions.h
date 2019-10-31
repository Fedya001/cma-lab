#pragma once

template<class T>
struct LDLTDecomposition {
  std::vector<bool> diagonal;
  SquareMatrix<T> low;
};

template<class T>
struct DLUDecomposition {
  std::vector<size_t> rows_permutations;
  SquareMatrix<T> low_up;
};