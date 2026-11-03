cmake_minimum_required(VERSION 3.15...3.31)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

project(prodigy_ginac LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

# The wheel build (pyproject.toml) turns both off, so it needs neither googletest nor the submodules.
option(PRODIGY_GINAC_BUILD_TESTS "Build the C++ tests (needs the googletest submodule)" ON)
option(PRODIGY_GINAC_BUILD_DEMO "Build the C++ demo executable" ON)
# Directory with the license files of GiNaC, CLN and GMP; the wheel build installs them into the wheel.
set(PRODIGY_GINAC_THIRD_PARTY_LICENSES "" CACHE PATH "Directory with license files of bundled libraries")

# find ginac on system and add headers to include path
find_package(GiNaC REQUIRED)
include_directories(${GINAC_INCLUDE_DIRS})

# so that tests and bindings can easily include the library headers
include_directories(src)

# actual library
add_subdirectory(src)

# tests
if (PRODIGY_GINAC_BUILD_TESTS)
    add_subdirectory(lib/googletest)
    add_subdirectory(test)
endif()

# python bindings: pybind11 from the Python environment (pip install pybind11), else the submodule
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
find_package(pybind11 CONFIG QUIET)
if (NOT pybind11_FOUND)
    add_subdirectory(lib/pybind11)
endif()
pybind11_add_module(prodigy_ginac python/bindings.cpp)
target_link_libraries(prodigy_ginac PRIVATE ${CMAKE_PROJECT_NAME}_lib)
if (DEFINED SKBUILD_PROJECT_VERSION)
    target_compile_definitions(prodigy_ginac PRIVATE VERSION_INFO=${SKBUILD_PROJECT_VERSION})
endif()
if (CMAKE_SYSTEM_NAME STREQUAL "Linux")
    # keep the symbols of statically linked GiNaC, CLN and GMP out of the module's dynamic symbol table
    target_link_options(prodigy_ginac PRIVATE "LINKER:--exclude-libs,ALL")
endif()

install(TARGETS prodigy_ginac LIBRARY DESTINATION .)
if (PRODIGY_GINAC_THIRD_PARTY_LICENSES AND DEFINED SKBUILD_METADATA_DIR)
    install(DIRECTORY "${PRODIGY_GINAC_THIRD_PARTY_LICENSES}/"
            DESTINATION "${SKBUILD_METADATA_DIR}/licenses/third_party")
endif()
