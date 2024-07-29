#!/bin/bash

mkdir -p $BUILD_PREFIX/src/cpp_module
cp -r $RECIPE_DIR/../../src/cpp_module/* $BUILD_PREFIX/src/cpp_module/

$PYTHON setup.py install
