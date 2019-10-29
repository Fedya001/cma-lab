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

  void SwapRows(size_t first_row, size_t second_row);
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

template<class M>
void PrintMatrix(std::ostream& out, const SquareMatrix<M>& matrix,
                 bool latex_syntax_on = false, uint8_t precision = 2,
                 const std::string& label = "A") {
  out << std::fixed << std::setprecision(precision);
  out << label << " = " <<
      (latex_syntax_on ? "\\left(\\begin{array}{ *{" + std::to_string(matrix.GetDim()) + "}{c} }\n" : "{\n");

  bool first_row = true;
  for (const auto& row : matrix.data_) {
    if (!first_row) {
      out << (latex_syntax_on ? " \\\\\n" : "\n");
    }

    bool first_in_a_row = true;
    for (const auto& element : row) {
      if (!first_in_a_row) {
        out << (latex_syntax_on ? " & " : "\t");
      }
      out << element;
      first_in_a_row = false;
    }
    first_row = false;
  }
  out << (latex_syntax_on ? "\n\\end{array}\\right)" : "\n}");
}
