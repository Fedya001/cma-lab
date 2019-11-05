#pragma once

#include "square_matrix.h"

#include <experimental/filesystem>
#include <fstream>
#include <iostream>

namespace loader {

namespace filesystem = std::experimental::filesystem;

inline void CheckFileExistence(const std::string& filename_path) {
  if (!filesystem::exists(filename_path) || !filesystem::is_regular_file(filename_path)) {
    throw std::logic_error("'filename_path' doesn't exist or is not a file");
  }
}

inline void CheckDirectoryExistence(const std::string& directory_path) {
  if (!filesystem::exists(directory_path) || !filesystem::is_directory(directory_path)) {
    throw std::logic_error("'directory_path' doesn't exists or is not a directory");
  }
}

std::vector<std::string> ListFiles(const std::string& directory_path);

template<class T>
SquareMatrix<T> LoadMatrix(std::istream& input) {
  size_t dim;
  input >> dim;

  T element;
  SquareMatrix<T> matrix(dim);
  for (size_t row = 0; row < dim; ++row) {
    for (size_t column = 0; column < dim; ++column) {
      input >> element;
      matrix[row][column] = element;
    }
  }

  return matrix;
}

template<class T>
std::vector<T> LoadVector(std::istream& input) {
  size_t dim;
  input >> dim;

  std::vector<T> vector(dim);
  for (auto& element : vector) {
    input >> element;
  }
  return vector;
}

template<class T>
SquareMatrix<T> LoadMatrix(const std::string& filename_path) {
  CheckFileExistence(filename_path);
  std::ifstream input(filename_path);
  return LoadMatrix<T>(input);
}

template<class T>
std::vector<T> LoadVector(const std::string& filename_path) {
  CheckFileExistence(filename_path);
  std::ifstream input(filename_path);
  return LoadVector<T>(input);
}

template<class T>
System<T> LoadSystem(const std::string& filename_path) {
  CheckFileExistence(filename_path);
  std::ifstream input(filename_path);
  return {LoadMatrix<T>(input), LoadVector<T>(input)};
}

template<class T>
std::vector<SquareMatrix<T>> LoadMatrices(const std::string& directory_path) {
  auto files = ListFiles(directory_path);

  std::vector<SquareMatrix<T>> matrices;
  matrices.reserve(files.size());
  for (const auto& file : files) {
    matrices.push_back(LoadMatrix<T>(file));
  }
  return matrices;
}

template<class T>
std::vector<std::vector<T>> LoadVectors(const std::string& directory_path) {
  auto files = ListFiles(directory_path);

  std::vector<std::vector<T>> vectors;
  vectors.reserve(files.size());
  for (const auto& file : files) {
    vectors.push_back(LoadVector<T>(file));
  }
  return vectors;
}

template<class T>
std::vector<System<T>> LoadSystems(const std::string& directory_path) {
  auto files = ListFiles(directory_path);

  std::vector<System<T>> systems;
  systems.reserve(files.size());
  for (const auto& file : files) {
    systems.push_back(LoadSystem<T>(file));
  }
  return systems;
}

} // namespace loader
