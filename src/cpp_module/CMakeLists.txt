cmake_minimum_required(VERSION 3.14)  # 设置最低的CMake版本
project(alexander_poly)                   # 定义项目名称

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake_file")

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
set(CMAKE_LINK_SEARCH_END_STATIC FALSE)
set(PYBIND11_FINDPYTHON OFF)

# 寻找Python环境
set (Python_FIND_VIRTUALENV ONLY)
set(CMAKE_FIND_ROOT_PATH_MODE_PACKAGE NEVER)
# find_package(Python COMPONENTS Development NumPy REQUIRED)
# message(STATUS "Python_EXECUTABLE: ${Python_EXECUTABLE}")

# 寻找pybind11库

find_package(pybind11 CONFIG REQUIRED)

# Set eigen3 path
set(EIGEN3_INCLUDE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/Eigen/")

# find_package(GMP REQUIRED)
# message(STATUS "GMP Include Directory: ${GMP_LIBRARIES}")

# find package cln and ginac
set (CLN_DIR "${CMAKE_CURRENT_SOURCE_DIR}/cln/build/")
find_package(CLN REQUIRED)

# if(CLN_FOUND)
#   message(STATUS "Found CLN: Yes")
#   message(STATUS "CLN Include Directory: ${CLN_INCLUDE_DIRS}")
#   message(STATUS "CLN Libraries: ${CLN_LIBRARIES}")
#   message(STATUS "CLN Version: ${CLN_VERSION}")
# else()
#   message(STATUS "Found CLN: No")
# endif()

# set 优先找so文件
# find_package(GINAC REQUIRED)

set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a" ".lib")
# message(${CMAKE_CURRENT_SOURCE_DIR})

find_library(GINAC_LIBRARY 
  NAMES ginac
  PATHS "${CMAKE_CURRENT_SOURCE_DIR}/ginac/build/ginac/"
  NO_DEFAULT_PATH
)
# find_library(CLN_LIBRARY
#   NAMES cln
#   PATHS "${CMAKE_CURRENT_SOURCE_DIR}/cln/build/lib/"
#   NO_DEFAULT_PATH
# )

if(NOT GINAC_LIBRARY)
  message(FATAL_ERROR "GINAC library not found")
else()
  message(STATUS "Found GINAC: Yes")
  message(STATUS "GINAC Library: ${GINAC_LIBRARY}")
endif()

# if(NOT CLN_LIBRARY)
#   message(FATAL_ERROR "CLN library not found")
# else()
#   message(STATUS "Found CLN: Yes")
#   message(STATUS "CLN Library: ${CLN_LIBRARY}")
# endif()

file(GLOB SOURCE_FILES "src/*.cpp")
# file(GLOB HEADER_FILES "src/*.h")
# 创建一个模块
pybind11_add_module(alexander_poly MODULE ${SOURCE_FILES})

target_link_libraries(alexander_poly PRIVATE ${GINAC_LIBRARY} ${CLN_LIBRARIES})

include_directories(${GMP_INCLUDE_DIRS} ${EIGEN3_INCLUDE_DIR})
# install(TARGETS alexander_poly
#   LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
#   RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
# )

# install(FILES ${GINAC_LIBRARY}
#   DESTINATION ${CMAKE_INSTALL_LIBDIR}
# )