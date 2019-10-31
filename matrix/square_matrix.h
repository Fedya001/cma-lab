#pragma once

#include <iomanip>
#include <sstream>
#include <vector>

template<class T>
class SquareMatrix {
 public:
  explicit SquareMatrix(size_t dim);

  explicit SquareMatrix(const std::vector<std::vector<T>>& data);
  explicit SquareMatrix(std::vector<std::vector<T>>&& data);

  SquareMatrix(const SquareMatrix& other) = default;
  SquareMatrix(SquareMatrix&& other) noexcept = default;

  SquareMatrix<T>& operator=(const SquareMatrix<T>& other);
  SquareMatrix<T>& operator=(SquareMatrix<T>&& other) noexcept;

  SquareMatrix<T>& operator*=(const SquareMatrix<T>& other);

  [[nodiscard]] std::vector<T>& operator[](size_t index);
  [[nodiscard]] const std::vector<double>& at(size_t index) const;

  void MultiplyRow(size_t row_index, T coefficient);
  void MultiplyColumn(size_t row_index, T coefficient);

  void SwapRows(size_t first_row, size_t second_row, size_t suffix = 0);
  void SwapColumns(size_t first_column, size_t second_column, size_t suffix = 0);

  [[nodiscard]] size_t GetDim() const;

  template<class M>
  friend SquareMatrix<M> operator*(const SquareMatrix<M>& lhs, const SquareMatrix<M>& rhs);
  template<class M>
  friend void PrintMatrix(std::ostream& out, const SquareMatrix<M>& matrix,
                          bool latex_syntax_on, uint8_t precision, const std::string& label);

 protected:
  size_t dim_;
  std::vector<std::vector<T>> data_;

 private:
  void EnsureMatrixStructure();
};

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
void SquareMatrix<T>::SwapRows(size_t first_row, size_t second_row, size_t suffix) {
  if (suffix == 0) {
    std::swap(data_[first_row], data_[second_row]);
  } else {
    for (size_t column_index = suffix; column_index < dim_; ++column_index) {
      std::swap(data_[first_row][column_index], data_[second_row][column_index]);
    }
  }
}

template<class T>
void SquareMatrix<T>::SwapColumns(size_t first_column, size_t second_column, size_t suffix) {
  std::vector<T> buffer;
  buffer.reserve(dim_);

  for (size_t row_index = suffix; row_index < dim_; ++row_index) {
    buffer.push_back(data_[row_index][first_column]);
    data_[row_index][first_column] = data_[row_index][second_column];
  }

  for (size_t row_index = suffix; row_index < dim_; ++row_index) {
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

template<class M>
SquareMatrix<M> operator*(const SquareMatrix<M>& lhs, const SquareMatrix<M>& rhs) {
  if (lhs.GetDim() != rhs.GetDim()) {
    throw std::logic_error("Can't multiply incompatible matrices");
  }
  std::vector<std::vector<M>> data;
  data.reserve(lhs.dim_);

  std::vector<M> line;
  for (size_t row = 0; row < lhs.dim_; ++row) {
    line.clear();
    for (size_t column = 0; column < rhs.dim_; ++column) {
      M sum = M();
      for (size_t index = 0; index < lhs.dim_; ++index) {
        sum += lhs.at(row).at(index) * rhs.at(index).at(column);
      }
      line.push_back(sum);
    }
    data.push_back(line);
  }

  return SquareMatrix<M>(std::move(data));
}
