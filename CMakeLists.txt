# -----------------------------------------------------------------------------
#
# Copyright (C) 2021 CERN & University of Surrey for the benefit of the
# BioDynaMo collaboration. All Rights Reserved.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
#
# See the LICENSE file distributed with this work for details.
# See the NOTICE file distributed with this work for additional information
# regarding copyright ownership.
#
# -----------------------------------------------------------------------------
cmake_minimum_required(VERSION 3.19.3)
project(pancreatic_tumor_new)

set(CMAKE_CXX_STANDARD 17)

find_package(BioDynaMo REQUIRED)
include(${BDM_USE_FILE})

# Add src/include to the public include path so you can do: #include "ptm/..."
# (You asked to keep headers under src/include, not top-level include/)
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/include")
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/src/")

# Collect sources/headers
file(GLOB_RECURSE PROJECT_HEADERS
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.h"
     "${CMAKE_CURRENT_SOURCE_DIR}/src/include/**/*.h")
file(GLOB_RECURSE PROJECT_SOURCES
     "${CMAKE_CURRENT_SOURCE_DIR}/src/*.cc")

# Build executable the BioDynaMo way (as you already do)
bdm_add_executable(${CMAKE_PROJECT_NAME}
                   HEADERS ${PROJECT_HEADERS}
                   SOURCES ${PROJECT_SOURCES}
                   LIBRARIES ${BDM_REQUIRED_LIBRARIES})

# Consider all files in test/ for GoogleTests.
# include_directories("test")
# file(GLOB_RECURSE TEST_SOURCES test/*.cc)
# file(GLOB_RECURSE TEST_HEADERS test/*.h)

# bdm_add_test(${CMAKE_PROJECT_NAME}-test
#              SOURCES ${TEST_SOURCES}
#              HEADERS ${TEST_HEADERS}
#              LIBRARIES ${BDM_REQUIRED_LIBRARIES} ${CMAKE_PROJECT_NAME})
