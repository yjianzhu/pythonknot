cmake_minimum_required(VERSION 3.12)
project(pdb_parser)

# 设置C++标准
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# 查找Python
set (PYBIND11_FINDPYTHON OFF)
set (Python_FIND_VIRTUALENV ONLY)
#find_package(Python COMPONENTS Interpreter Development REQUIRED)

# 查找pybind11
find_package(pybind11 REQUIRED)

# 添加源文件
set(SOURCES
    src/pdb.cpp
)

# 添加头文件目录
include_directories(
    ${CMAKE_CURRENT_SOURCE_DIR}/src
)

# 创建Python模块
pybind11_add_module(pdb_parser ${SOURCES})

# 添加编译选项
target_compile_options(pdb_parser PRIVATE
    -Wall
    -O3
)

# 安装配置（可选）
# install(TARGETS pdb_parser DESTINATION ${CMAKE_INSTALL_PREFIX})