cmake_minimum_required(VERSION 3.4...3.18)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

project(pygin)

set(CMAKE_CXX_STANDARD 14)

# find ginac on system and add headers to include path
find_package(GiNaC REQUIRED)
include_directories(${GINAC_INCLUDE_DIRS})

# so that tests and bindings can easily include the library headers
include_directories(src)

# actual library
add_subdirectory(src)

# tests
add_subdirectory(lib/googletest)
add_subdirectory(test)

# python bindings
add_subdirectory(lib/pybind11)
pybind11_add_module(pygin python/bindings.cpp)
target_link_libraries(pygin PUBLIC ${CMAKE_PROJECT_NAME}_lib)
target_compile_definitions(pygin PRIVATE VERSION_INFO=${EXAMPLE_VERSION_INFO})
