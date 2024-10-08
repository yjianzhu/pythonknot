cmake_minimum_required(VERSION 3.14)  # 设置最低的CMake版本
project(alexander_poly)                   # 定义项目名称

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_file")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_LINK_SEARCH_END_STATIC FALSE)
set(PYBIND11_FINDPYTHON OFF)

find_package(GMP REQUIRED)
message(STATUS "GMP Include Directory: ${GMP_LIBRARIES}")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".a" ".so" ".lib")

# 寻找Python环境
set (Python_FIND_VIRTUALENV ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE NEVER)
# find_package(Python COMPONENTS Development NumPy REQUIRED)
# message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")

# 寻找pybind11库

find_package(pybind11 CONFIG REQUIRED)

# Set eigen3 path
# Find Eigen3, if not found, set the path manually
find_package (Eigen3 3.3 NO_MODULE)
if (NOT Eigen3_FOUND)
  message(STATUS "Eigen3 not found, set the path manually")
  set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/eigen3/")
  message(STATUS "Eigen3 Include Directory: ${EIGEN3_INCLUDE_DIR}")
else()
  message(STATUS "Eigen3 Include Directory: ${EIGEN3_INCLUDE_DIR}")
endif()


# find package cln and ginac
set (CLN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cln/build/")
find_package(CLN REQUIRED)

if(CLN_FOUND)
  message(STATUS "Found CLN: Yes")
  message(STATUS "CLN Include Directory: ${CLN_INCLUDE_DIR}")
  message(STATUS "CLN Libraries: ${CLN_LIBRARIES}")
  message(STATUS "CLN Version: ${CLN_VERSION}")
else()
  message(STATUS "Found CLN: No")
endif()

# set 优先找so文件
# find_package(GINAC REQUIRED)

# message(${CMAKE_CURRENT_SOURCE_DIR})

find_library(GINAC_LIBRARY 
  NAMES ginac
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/ginac/install/lib"
  NO_DEFAULT_PATH
)

find_path(GINAC_INCLUDE_DIR
  NAMES ginac/ginac.h
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/ginac/install/include"
  NO_DEFAULT_PATH
)

if(NOT GINAC_LIBRARY)
  message(FATAL_ERROR "GINAC library not found")
else()
  message(STATUS "Found GINAC: Yes")
  message(STATUS "GINAC Library: ${GINAC_LIBRARY}")
  message(STATUS "GINAC Include Directory: ${GINAC_INCLUDE_DIR}")
endif()

file(GLOB SOURCE_FILES "src/*.cpp")
# file(GLOB HEADER_FILES "src/*.h")
# 创建一个模块
pybind11_add_module(alexander_poly MODULE ${SOURCE_FILES})

target_link_libraries(alexander_poly PRIVATE ${GINAC_LIBRARY} ${CLN_LIBRARIES} ${GMP_LIBRARIES})

# if find Eigen3, include the directory
if (Eigen3_FOUND)
  target_link_libraries(alexander_poly PRIVATE Eigen3::Eigen)
else()
  message(STATUS "Include Eigen3 manually")
  include_directories(${EIGEN3_INCLUDE_DIR})
endif()

include_directories(${GINAC_INCLUDE_DIR} ${CLN_INCLUDE_DIR})