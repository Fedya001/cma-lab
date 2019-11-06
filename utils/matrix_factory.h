#pragma once

#include <chrono>
#include <random>

template<class T>
class MatrixFactory {
 private:
  template<bool is_integral, class M>
  struct uniform_distribution_selector;

  template<class M>
  struct uniform_distribution_selector<true, M> {
    using type = typename std::uniform_int_distribution<T>;
  };

  template<class M>
  struct uniform_distribution_selector<false, M> {
    using type = typename std::uniform_real_distribution<T>;
  };

 public:
  using uniform_distribution_type = typename uniform_distribution_selector<std::is_integral<T>::value, T>::type;

  explicit MatrixFactory(T lower_bound = std::numeric_limits<T>::min(),
                         T upper_bound = std::numeric_limits<T>::max(),
                         long long seed = std::chrono::duration_cast<std::chrono::microseconds>
                             (std::chrono::system_clock::now().time_since_epoch()).count());

  SquareMatrix<T> CreateRandomMatrix(size_t dim, bool symmetric = false);
  SquareMatrix<T> CreateLowdiagMatrix(size_t dim);

  std::vector<T> CreateRandomVector(size_t length);

  static SquareMatrix<T> CreateIdentity(size_t dim);

 private:
  std::mt19937 generator_;
  uniform_distribution_type distribution_;
};

template<class T>
MatrixFactory<T>::MatrixFactory(T lower_bound, T upper_bound, long long int seed)
    : generator_(seed),
      distribution_(lower_bound, upper_bound) {}

template<class T>
SquareMatrix<T> MatrixFactory<T>::CreateRandomMatrix(size_t dim, bool symmetric) {
  std::vector<std::vector<T>> data;
  data.reserve(dim);

  for (size_t row_index = 0; row_index < dim; ++row_index) {
    std::vector<T> row;
    row.reserve(dim);

    if (symmetric) {
      row.resize(row_index, T(0));
    }

    for (size_t column_index = (symmetric ? row_index : 0); column_index < dim; ++column_index) {
      row.push_back(distribution_(generator_));
    }

    if (symmetric) {
      for (size_t column_index = 0; column_index < row_index; ++column_index) {
        row[column_index] = data[column_index][row_index];
      }
    }

    data.push_back(row);
  }

  return SquareMatrix<T>(std::move(data));
}

template<class T>
SquareMatrix<T> MatrixFactory<T>::CreateLowdiagMatrix(size_t dim) {
  std::vector<std::vector<T>> data;
  data.reserve(dim);

  for (size_t row_index = 0; row_index < dim; ++row_index) {
    std::vector<T> row;
    size_t length = std::min(dim, row_index + 2);
    row.reserve(length);
    for (size_t column_index = 0; column_index < length; ++column_index) {
      row.push_back(distribution_(generator_));
    }
    row.resize(dim, T(0));
    data.push_back(row);
  }

  return SquareMatrix<T>(std::move(data));
}

template<class T>
std::vector<T> MatrixFactory<T>::CreateRandomVector(size_t length) {
  std::vector<T> result;
  for (size_t index = 0; index < length; ++index) {
    result.push_back(distribution_(generator_));
  }
  return result;
}

template<class T>
SquareMatrix<T> MatrixFactory<T>::CreateIdentity(size_t dim) {
  std::vector<T> row(dim, T(0));
  std::vector<std::vector<T>> data;

  for (size_t row_index = 0; row_index < dim; ++row_index) {
    row[row_index] = T(1);
    data.push_back(row);
    row[row_index] = T(0);
  }

  return SquareMatrix<T>(std::move(data));
}
