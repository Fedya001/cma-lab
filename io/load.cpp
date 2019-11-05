#include "load.h"

namespace loader {

std::vector<std::string> ListFiles(const std::string& directory_path) {
  CheckDirectoryExistence(directory_path);

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

  return files;
}
} // namespace loader