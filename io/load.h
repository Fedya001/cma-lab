#pragma once

#include "square_matrix.h"

#include <experimental/filesystem>
#include <fstream>
#include <iostream>

namespace loader {

namespace filesystem = std::experimental::filesystem;

template<class T>
SquareMatrix<T> LoadMatrix(const std::string& filename_path) {
  if (!filesystem::exists(filename_path)) {
    throw std::invalid_argument("LoadMatrix(): 'filename_path' doesn't exists");
  }
  if (!filesystem::is_regular_file(filename_path)) {
    throw std::invalid_argument("LoadMatrix(): 'filename_path' is not a file");
  }

  std::ifstream input(filename_path);
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
std::vector<SquareMatrix<T>> LoadMatrices(const std::string& directory_path) {
  if (!filesystem::exists(directory_path)) {
    throw std::invalid_argument("LoadMatrices(): 'directory_path' doesn't exists");
  }
  if (!filesystem::is_directory(directory_path)) {
    throw std::invalid_argument("LoadMatrices(): 'directory_path' is not a directory");
  }

  std::vector<std::string> files;
  try {
    filesystem::recursive_directory_iterator iter(directory_path);
    // Iterate till end
    while (iter != filesystem::recursive_directory_iterator()) {
      if (filesystem::is_directory(iter->path())) {
        // c++17 Filesystem API to skip current directory iteration
        iter.disable_recursion_pending();
      } else {
        // Add the name in vector
        files.push_back(iter->path().string());
      }

      std::error_code error;
      iter.increment(error);
      if (error) {
        std::cerr << "Error While Accessing : " << iter->path().string()
                  << " :: " << error.message() << '\n';
      }
    }
  }
  catch (std::system_error& e) {
    std::cerr << e.what();
  }

  std::vector<SquareMatrix<T>> matrices;
  matrices.reserve(files.size());

  for (const auto& file : files) {
    matrices.push_back(LoadMatrix<T>(file));
  }

  return matrices;
}

} // namespace loader
