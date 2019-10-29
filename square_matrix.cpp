#include "square_matrix.h"

template<class T>
SquareMatrix<T>::SquareMatrix(size_t dim)
    : dim_(dim),
      data_(dim, std::vector<T>(dim)) {}

template<class T>
SquareMatrix<T>::SquareMatrix(const std::vector<std::vector<T>>& data)
    : dim_(data.size()),
      data_(data) {
  EnsureMatrixStructure();
}

template<class T>
SquareMatrix<T>::SquareMatrix(std::vector<std::vector<T>>&& data)
    : dim_(data.size()),
      data_(std::move(data)) {
  data = std::vector<std::vector<T>>();
  EnsureMatrixStructure();
}

template<class T>
SquareMatrix<T>& SquareMatrix<T>::operator=(const SquareMatrix<T>& other) {
  dim_ = other.dim_;
  data_ = other.data_;
  return *this;
}

template<class T>
SquareMatrix<T>& SquareMatrix<T>::operator=(SquareMatrix<T>&& other) noexcept {
  dim_ = other.dim_;
  data_ = std::move(other.data_);
  return *this;
}

template<class T>
SquareMatrix<T>& SquareMatrix<T>::operator*=(const SquareMatrix<T>& other) {
  return *this = *this * other;
}

template<class T>
std::vector<T>& SquareMatrix<T>::operator[](size_t index) {
  if (index >= data_.size()) {
    throw std::out_of_range("Invalid index in operator[]");
  }
  return data_[index];
}

template<class T>
const std::vector<double>& SquareMatrix<T>::at(size_t index) const {
  if (index >= data_.size()) {
    throw std::out_of_range("Invalid index in operator[]");
  }
  return data_.at(index);
}

template<class T>
void SquareMatrix<T>::MultiplyRow(size_t row_index, T coefficient) {
  if (row_index >= data_.size()) {
    throw std::out_of_range("Invalid row index");
  }
  auto& row = data_[row_index];
  for (T& element : row) {
    element *= coefficient;
  }
}

template<class T>
void SquareMatrix<T>::MultiplyColumn(size_t column_index, T coefficient) {
  if (column_index >= data_.size()) {
    throw std::out_of_range("Invalid column index");
  }
  for (size_t row_index = 0; row_index < dim_; ++row_index) {
    data_[row_index][column_index] *= coefficient;
  }
}

template<class T>
void SquareMatrix<T>::SwapRows(size_t first_row, size_t second_row) {
  std::swap(data_[first_row], data_[second_row]);
}

template<class T>
void SquareMatrix<T>::SwapColumns(size_t first_column, size_t second_column) {
  std::vector<T> buffer;
  buffer.reserve(dim_);

  for (size_t row_index = 0; row_index < dim_; ++row_index) {
    buffer.push_back(data_[row_index][first_column]);
    data_[row_index][first_column] = data_[row_index][second_column];
  }

  for (size_t row_index = 0; row_index < dim_; ++row_index) {
    data_[row_index][second_column] = buffer[row_index];
  }
}

template<class T>
size_t SquareMatrix<T>::GetDim() const {
  return dim_;
}

template<class T>
void SquareMatrix<T>::EnsureMatrixStructure() {
  for (const auto& row : data_) {
    if (row.size() != dim_) {
      throw std::invalid_argument("Invalid matrix structure");
    }
  }
}

// explicit instantiation of template
template
class SquareMatrix<double>;