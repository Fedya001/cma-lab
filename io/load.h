#pragma once

#include "square_matrix.h"
#include "three_diagonal_matrix.h"

#include <experimental/filesystem>
#include <fstream>
#include <iostream>

// Because the lack of time, now there is some code duplication for each
// type of matrix and system
// TODO:
//  Expand logic for inputting different types of matrices and systems
//  instead of simple duplicating code for each type

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
ThreeDiagonalMatrix<T> LoadThreeDiagonalMatrix(std::istream& input) {
  int32_t dim;
  input >> dim;
  if (dim <= 0) {
    throw std::logic_error("Matrix dimension must be positive");
  }

  if (dim == 1) {
    T element;
    input >> element;
    return ThreeDiagonalMatrix<T>({{0, element, 0}});
  } else {
    std::vector<std::vector<T>> data;
    data.reserve(dim);

    std::vector<T> boundary_row(3);
    input >> boundary_row[1] >> boundary_row[2];
    data.push_back(boundary_row);

    for (int32_t index = 0; index < dim - 2; ++index) {
      std::vector<T> row(3);
      input >> row[0] >> row[1] >> row[2];
      data.push_back(row);
    }

    input >> boundary_row[0] >> boundary_row[1];
    boundary_row[2] = T();
    data.push_back(boundary_row);

    return ThreeDiagonalMatrix(std::move(data));
  }
}

template<class T>
SquareMatrix<T> LoadMatrix(std::istream& input) {
  int32_t dim;
  input >> dim;
  if (dim <= 0) {
    throw std::logic_error("Matrix dimension must be positive");
  }

  T element;
  SquareMatrix<T> matrix(dim);
  for (int32_t row = 0; row < dim; ++row) {
    for (int32_t column = 0; column < dim; ++column) {
      input >> element;
      matrix[row][column] = element;
    }
  }

  return matrix;
}

template<class T>
std::vector<T> LoadVector(std::istream& input) {
  int32_t dim;
  input >> dim;
  if (dim <= 0) {
    throw std::logic_error("Vector dimension must be positive");
  }

  std::vector<T> vector(dim);
  for (auto& element : vector) {
    input >> element;
  }
  return vector;
}

template<class T>
ThreeDiagonalMatrix<T> LoadThreeDiagonalMatrix(const std::string& filename_path) {
  CheckFileExistence(filename_path);
  std::ifstream input(filename_path);
  return LoadThreeDiagonalMatrix<T>(input);
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
ThreeDiagonalSystem<T> LoadThreeDiagonalSystem(const std::string& filename_path) {
  CheckFileExistence(filename_path);
  std::ifstream input(filename_path);
  auto matrix = LoadThreeDiagonalMatrix<T>(input);
  auto column = LoadVector<T>(input);
  return {matrix, column};
}

template<class T>
System<T> LoadSystem(const std::string& filename_path) {
  CheckFileExistence(filename_path);
  std::ifstream input(filename_path);
  auto matrix = LoadMatrix<T>(input);
  auto column = LoadVector<T>(input);
  return {matrix, column};
}

template<class T>
std::vector<SquareMatrix<T>> LoadThreeDiagonalMatrices(const std::string& directory_path) {
  auto files = ListFiles(directory_path);

  std::vector<ThreeDiagonalMatrix<T>> matrices;
  matrices.reserve(files.size());
  for (const auto& file : files) {
    matrices.push_back(LoadThreeDiagonalMatrices<T>(file));
  }
  return matrices;
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
std::vector<ThreeDiagonalSystem<T>> LoadThreeDiagonalSystems(const std::string& directory_path) {
  auto files = ListFiles(directory_path);

  std::vector<ThreeDiagonalSystem<T>> systems;
  systems.reserve(files.size());
  for (const auto& file : files) {
    systems.push_back(LoadThreeDiagonalSystem<T>(file));
  }
  return systems;
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
