cmake_minimum_required(VERSION 3.10)
project(cma-lab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

include_directories(matrix utils)
add_executable(cma-lab main.cpp)
target_link_libraries(cma-lab stdc++fs)