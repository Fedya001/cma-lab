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
T ComputeDiffEuclidNorm(int32_t size, const std::tuple<T, T, T>& lhs, const std::tuple<T, T, T>& rhs) {
  return std::sqrt(
      std::pow((std::get<0>(lhs) - std::get<0>(rhs)), 2.0) +
      std::pow((std::get<1>(lhs) - std::get<1>(rhs)), 2.0) * (size - 2) +
      std::pow((std::get<2>(lhs) - std::get<2>(rhs)), 2.0)
  );
}

template<class T>
std::tuple<T, T, T> SolveSystem(int32_t size, T omega, T accuracy = T(1e-9)) {
  if (size <= 0) {
    throw std::logic_error("Size of the system must be positive");
  }

  if (size == 1) {
    return {T(1), T(1), T(1)};
  }

  std::tuple<T, T, T> current{0, 0, 0}, next{0, 0, 0};

  T one_minus_omega = T(1) - omega;
  do {
    std::get<0>(next) = one_minus_omega * std::get<0>(current) + omega *
        (T(1) - (size - 2) * std::get<1>(current) - std::get<2>(current)) / size;
    std::get<1>(next) = one_minus_omega * std::get<1>(current) + omega *
        (T(1) - std::get<0>(next) - std::get<2>(current)) / size;
    std::get<2>(next) = one_minus_omega * std::get<2>(current) + omega *
        (T(1) - std::get<0>(next) - (size - 2) * std::get<1>(next)) / size;
    current.swap(next);
  } while (ComputeDiffEuclidNorm(size, current, next) > accuracy);
  return current;
}

