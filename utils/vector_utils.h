#pragma once

#include <cmath>
#include <stdexcept>

template<class T>
std::vector<T> NoiseVector(const std::vector<T>& vector, T noise_radius = 1.0) {
  std::mt19937_64 gen(std::random_device{}());
  std::uniform_real_distribution<T> distr(-noise_radius, noise_radius);

  auto copy = vector;
  for (T& element : copy) {
    element += distr(gen);
  }
  return copy;
}

template<class T>
T GetDiffEuclidNorm(const std::vector<T>& lhs, const std::vector<T>& rhs) {
  if (lhs.size() != rhs.size()) {
    throw std::invalid_argument("GetDiffEuclidNorm(): vector sized must be equal");
  }

  T sum = T();
  for (int32_t index = 0; index < static_cast<int32_t>(lhs.size()); ++index) {
    T diff = lhs[index] - rhs[index];
    sum += diff * diff;
  }
  return std::sqrt(sum);
}