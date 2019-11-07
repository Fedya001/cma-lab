#pragma once

#include <cmath>
#include <vector>
#include <stdexcept>

// TODO: Make a template matrix hierarchy

template<class T>
class ThreeDiagonalMatrix {
 public:
  explicit ThreeDiagonalMatrix(const std::vector<std::vector<T>>& data);
  explicit ThreeDiagonalMatrix(std::vector<std::vector<T>>&& data);

  std::vector<T> SolveSystem(std::vector<T> result, bool swap_rows = false) const;

  [[nodiscard]] std::vector<T>& operator[](int32_t index);
  [[nodiscard]] const std::vector<T>& operator[](int32_t index) const;

  [[nodiscard]] int32_t GetDim() const;
  [[nodiscard]] bool IsEmpty() const;

 private:
  int32_t dim_;
  // The first and the last row must have
  // trailing zeros
  std::vector<std::vector<T>> data_;

  void EnsureMatrixStructure() const;
  void OffsetVectorLeft(std::vector<T>& vector) const;
};

template<class T>
ThreeDiagonalMatrix<T>::ThreeDiagonalMatrix(const std::vector<std::vector<T>>& data)
    : dim_(data.size()),
      data_(data) {
  EnsureMatrixStructure();
}

template<class T>
ThreeDiagonalMatrix<T>::ThreeDiagonalMatrix(std::vector<std::vector<T>>&& data)
    : dim_(data.size()),
      data_(std::move(data)) {
  data = std::vector<std::vector<T>>();
  EnsureMatrixStructure();
}

template<class T>
std::vector<T> ThreeDiagonalMatrix<T>::SolveSystem(std::vector<T> result, bool swap_rows) const {
  // Make a copy to leave state unchanged
  auto data = data_;

  // Stage 1. Go down. Transform matrix to three-diagonal upper-triangular
  // with ones on the diagonal
  for (int32_t row = 0; row < dim_ - 1; ++row) {
    OffsetVectorLeft(data[row]);
    if (swap_rows && std::abs(data[row][0]) < std::abs(data[row + 1][0])) {
      data[row].swap(data[row + 1]);
      std::swap(result[row], result[row + 1]);
    }

    if (data[row].front() == 0) {
      std::string message = "Unable to solve system with sweep method.";
      if (!swap_rows) {
        message += " Consider set swap_rows to true.";
      } else {
        message += " Matrix is degenerate.";
      }
      throw std::logic_error(message);
    }

    result[row] /= data[row].front();
    for (int32_t index = 2; index >= 0; --index) {
      data[row][index] /= data[row].front();
    }

    result[row + 1] -= data[row + 1].front() * result[row];
    for (int32_t index = 2; index >= 0; --index) {
      data[row + 1][index] -= data[row + 1].front() * data[row][index];
    }
  }
  OffsetVectorLeft(data.back());
  result[dim_ - 1] /= data[dim_ - 1].front();
  data[dim_ - 1].front() = T(1);

  // Stage 2. Go up. Solve the system
  for (int32_t row = dim_ - 1; row >= 0; --row) {
    for (int32_t index : {1, 2}) {
      if (row >= index) {
        result[row - index] -= result[row] * data[row - index][index];
      }
    }
  }

  return result;
}

template<class T>
std::vector<T>& ThreeDiagonalMatrix<T>::operator[](int32_t index) {
  if (index < 0 || index >= static_cast<int32_t>(data_.size())) {
    throw std::out_of_range("Invalid index in operator[]");
  }
  return data_[index];
}

template<class T>
const std::vector<T>& ThreeDiagonalMatrix<T>::operator[](int32_t index) const {
  if (index < 0 || index >= static_cast<int32_t>(data_.size())) {
    throw std::out_of_range("Invalid index in operator[]");
  }
  return data_[index];
}

template<class T>
int32_t ThreeDiagonalMatrix<T>::GetDim() const {
  return dim_;
}

template<class T>
bool ThreeDiagonalMatrix<T>::IsEmpty() const {
  return dim_ == 0;
}

template<class T>
void ThreeDiagonalMatrix<T>::EnsureMatrixStructure() const {
  const std::string message = "Invalid three-diagonal matrix structure";
  for (const auto& row : data_) {
    if (row.size() != 3) {
      throw std::invalid_argument(message);
    }
  }
}

template<class T>
void ThreeDiagonalMatrix<T>::OffsetVectorLeft(std::vector<T>& vector) const {
  if (vector.empty()) {
    return;
  }

  T front = vector.front();
  for (int32_t index = 0; index < static_cast<int32_t>(vector.size()) - 1; ++index) {
    vector[index] = vector[index + 1];
  }
  vector.back() = front;
}

// TODO: create a normal class for systems
template<class T>
using ThreeDiagonalSystem = std::pair<ThreeDiagonalMatrix<T>, std::vector<T>>;
