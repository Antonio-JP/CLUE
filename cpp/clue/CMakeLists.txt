cmake_minimum_required(VERSION 3.19...3.28)

project(CLUE
        LANGUAGES CXX
        DESCRIPTION "Constrained LUmping of differential Equations in C++ - A faster implementation of linear model reduction.")

# Setting up the C++ standard used
set(CMAKE_CXX_STANDARD 20)

# Check whether `modulename` is correctly cloned in the `extern` directory.
macro(CHECK_SUBMODULE_PRESENT modulename)
  if(NOT EXISTS "${PROJECT_SOURCE_DIR}/extern/${modulename}/CMakeLists.txt")
    message(
      FATAL_ERROR
        "${modulename} submodule not cloned properly. \
        Please run `git submodule update --init --recursive` \
        from the main project directory")
  endif()
endmacro()

check_submodule_present(mqt-core)

# Including Boost into the linking and looking
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Including my own libraries
add_subdirectory(src)

## Code for Testing
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
        include(CTest) # We include CTest which is part of CMake

        # We check if this is the main file since we don't usually want users of the library to execute tests
        if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
                add_subdirectory(tests)
        endif ()
        enable_testing()
endif()

# If this is the main build, build possible applications
if(CMAKE_PROJECT_NAME STREQUAL PROJECT_NAME)
  add_subdirectory(apps)
endif()