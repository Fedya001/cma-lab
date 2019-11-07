#pragma once

#include <cmath>
#include <cstdint>
#include <numeric>
#include <stdexcept>
#include <vector>

// This function solves system of linear equations like this:
// n 1 1 ... 1 1 1   1
// 1 n 0 ... 0 0 1   1
// 1 0 n ... 0 0 1   1
// ............... = 1
// 1 0 0 ... n 0 1   1
// 1 0 0 ... 0 n 1   1
// 1 1 1 ... 1 1 n   1
// using relaxation method with omega parameter

template<class T>
T ComputeDiffEuclidNorm(const std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size()) {
    throw std::logic_error("Vector sizes must be equal");
  }

  T sum = T();
  for (int32_t index = 0; index < static_cast<int32_t>(lhs.size()); ++index) {
    auto diff = lhs[index] - rhs[index];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}

template<class T>
std::vector<T> SolveSystem(int32_t size, T omega, T accuracy = T(1e-9)) {
  if (size <= 0) {
    throw std::logic_error("Size of the system must be positive");
  }

  if (size == 1) {
    return {T(1)};
  }

  std::vector<T> current(size, T(0)), next(size);

  T one_minus_omega = T(1) - omega;
  do {
    next[0] = one_minus_omega * current[0] + omega *
        (T(1) - std::accumulate(current.begin() + 1, current.end(), T())) / size;
    for (int32_t index = 1; index < size - 1; ++index) {
      next[index] = one_minus_omega * current[index] + omega *
          (T(1) - next[0] - current[size - 1]) / size;
    }
    next[size - 1] = one_minus_omega * current[size - 1] + omega *
        (T(1) - std::accumulate(next.begin(), next.end() - 1, T())) / size;

    current.swap(next);
  } while (ComputeDiffEuclidNorm(current, next) > accuracy);
  return current;
}

