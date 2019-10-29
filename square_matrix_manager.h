#pragma once

#include "square_matrix.h"

template<class T>
class SquareMatrixManager {
 public:
  enum class MatrixType {
    UPPER_DIAGONAL,
    LOWER_DIAGONAL,
    DIAGONAL,
    ORDINARY
  };

  // inner classes representing decompositions
  struct LDLTDecomposition {
    std::vector<bool> diagonal;
    SquareMatrix<T> low;
  };

  struct DLUDecomposition {
    std::vector<size_t> rows_permutations;
    SquareMatrix<T> low_up;
  };

  explicit SquareMatrixManager(SquareMatrix<T> matrix);

  SquareMatrix<T> GetMatrix() const;
  void SetMatrix(const SquareMatrix<T>& matrix);

  LDLTDecomposition PerformLDLT() const;
  DLUDecomposition PerformDLU(bool swap_rows = false) const;
  std::vector<T> SolveSystem(const std::vector<T>& result) const;

 private:
  SquareMatrix<T> matrix_;
};
