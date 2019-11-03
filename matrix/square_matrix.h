#pragma once

#include <iomanip>
#include <sstream>
#include <vector>
#include <unordered_set>
#include <unordered_map>

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
  [[nodiscard]] const std::vector<T>& operator[](size_t index) const;

  void MultiplyRow(size_t row_index, T coefficient);
  void MultiplyColumn(size_t row_index, T coefficient);

  void SwapRows(size_t first_row, size_t second_row);
  void SwapRows(std::vector<size_t> permutation);
  void SwapColumns(size_t first_column, size_t second_column);

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
const std::vector<T>& SquareMatrix<T>::operator[](size_t index) const {
  return data_[index];
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
  // std::vector.swap() takes O(1) time
  data_[first_row].swap(data_[second_row]);
}

template<class T>
void SquareMatrix<T>::SwapRows(std::vector<size_t> permutation) {
  // Empty permutation
  if (permutation.empty()) {
    return;
  }

  // Validate the permutation
  const std::string invalid_permutation_message = "Wrong permutation. ";
  if (permutation.size() != dim_) {
    throw std::invalid_argument(invalid_permutation_message +
        "Matrix dimension and permutation length don't match");
  }
  std::unordered_map<size_t, size_t> index_of(dim_);
  for (size_t index = 0; index < permutation.size(); ++index) {
    if (index_of.count(permutation[index])) {
      throw std::invalid_argument(invalid_permutation_message);
    }
    index_of[permutation[index]] = index;
  }

  if (index_of.size() != dim_) {
    throw std::invalid_argument(invalid_permutation_message);
  }

  // Smart in-place swaps
  // O(n) time, O(1) memory
  std::unordered_set<size_t> indices_left(permutation.begin(), permutation.end());
  size_t insert_index = 0;
  for (size_t operation = 0; operation < dim_; ++operation) {
    if (insert_index != permutation[insert_index]) {
      // Constant complexity swap
      data_[insert_index].swap(data_[permutation[insert_index]]);

      // Recover permutation invariance
      permutation[index_of.at(insert_index)] = permutation[insert_index];
      index_of[permutation[insert_index]] = index_of.at(insert_index);
      insert_index = permutation[insert_index];
      indices_left.erase(insert_index);
    } else {
      // Pull new index
      insert_index = *indices_left.begin();
      indices_left.erase(indices_left.begin());
    }
  }
}

template<class T>
void SquareMatrix<T>::SwapColumns(size_t first_column, size_t second_column) {
  std::vector<T> buffer;
  buffer.reserve(dim_);
  for (size_t row_index = 0; row_index <= dim_; ++row_index) {
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
        sum += lhs[row][index] * rhs[index][column];
      }
      line.push_back(sum);
    }
    data.push_back(line);
  }

  return SquareMatrix<M>(std::move(data));
}
