cmake_minimum_required(VERSION 3.14)  # 设置最低的CMake版本
project(homfly) 

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_file")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_LINK_SEARCH_END_STATIC FALSE)
set(PYBIND11_FINDPYTHON OFF)

# 寻找Python环境
set (Python_FIND_VIRTUALENV ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE NEVER)

find_package(pybind11 CONFIG REQUIRED)

set (GC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cpp_module/gc-8.2.8/install/")
find_package(LibGC REQUIRED)
if(LIBGC_FOUND)
  message(STATUS "Found LibGC: Yes")
  message(STATUS "LibGC Include Directory: ${LIBGC_INCLUDE_DIR}")
else()
  message(STATUS "Found LibGC: No")
endif()

message(STATUS "LibGC Include Directory: ${LIBGC_INCLUDE_DIR}")
message(STATUS "LibGC Library: ${LIBGC_LIBRARY}")

# Set eigen3 path
# Find Eigen3, if not found, set the path manually
find_package (Eigen3 3.3 NO_MODULE)
if (NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found, set the path manually")
  set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../cpp_module/eigen3/")
  message(STATUS "Eigen3 Include Directory: ${EIGEN3_INCLUDE_DIR}")
else()
  message(STATUS "Eigen3 Include Directory: ${EIGEN3_INCLUDE_DIR}")
endif()

file(GLOB SOURCE_FILES "src/*.cpp")
# file(GLOB HEADER_FILES "src/*.h")
# 创建一个模块
pybind11_add_module(homfly MODULE ${SOURCE_FILES})

# if find Eigen3, include the directory
if (Eigen3_FOUND)
  target_link_libraries(homfly PRIVATE Eigen3::Eigen)
else()
  message(STATUS "Include Eigen3 manually")
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

target_link_libraries(homfly PRIVATE ${LIBGC_LIBRARY})
target_include_directories(homfly PRIVATE ${LIBGC_INCLUDE_DIR})