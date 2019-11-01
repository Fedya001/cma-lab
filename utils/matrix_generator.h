#pragma once

#include <chrono>
#include <random>

template<class T>
class MatrixGenerator {
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

  explicit MatrixGenerator(T lower_bound = std::numeric_limits<T>::min(),
                           T upper_bound = std::numeric_limits<T>::max(),
                           const long long seed = std::chrono::duration_cast<std::chrono::microseconds>
                               (std::chrono::system_clock::now().time_since_epoch()).count())
      : generator_(seed),
        distribution_(lower_bound, upper_bound) {}

  SquareMatrix<T> generate(size_t dim) {
    std::vector<std::vector<T>> data;
    data.reserve(dim);

    for (size_t row_index = 0; row_index < dim; ++row_index) {
      std::vector<T> row;
      row.reserve(dim);
      for (size_t column_index = 0; column_index < dim; ++column_index) {
        row.push_back(distribution_(generator_));
      }
      data.push_back(row);
    }

    return SquareMatrix<T>(std::move(data));
  }

 private:
  std::mt19937 generator_;
  uniform_distribution_type distribution_;
};