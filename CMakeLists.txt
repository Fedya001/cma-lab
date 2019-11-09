cmake_minimum_required(VERSION 3.10)
project(cma-lab)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsanitize=address")

SET(CMAKE_CXX_FLAGS_DEBUG "-O3")
SET(CMAKE_CXX_FLAGS_RELEASE  "-O3")

include_directories(matrix utils io report)
add_executable(cma-lab main.cpp io/load.cpp report/report.cpp)
target_link_libraries(cma-lab stdc++fs)