cmake_minimum_required(VERSION 3.16)

project(CLUE)

set(CMAKE_CXX_STANDARD 20)

# Including Boost into the linking and looking
find_package(Boost)
include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

# Including my own libraries
include_directories("${CMAKE_CURRENT_SOURCE_DIR}/include")
include_directories(
        ${PROJECT_SOURCE_DIR}/include
        ${PROJECT_SOURCE_DIR}/src
)

file(GLOB all_SRCS
        "${PROJECT_SOURCE_DIR}/include/*.h"
        "${PROJECT_SOURCE_DIR}/include/*.hpp"
        "${PROJECT_SOURCE_DIR}/src/*.cpp"
        "${PROJECT_SOURCE_DIR}/src/*.c"
)

# Defining executable
add_executable(CLUE src/main.cpp src/RationalFunctions.cpp)
target_link_libraries(CLUE PUBLIC ${Boost_LIBRARIES})

## Code for Testing
include(CTest) # We include CTest which is part of CMake

if (${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR})
        add_subdirectory(tests)
endif ()

# We check if this is the main file
# you don't usually want users of your library to
# execute tests as part of their build

enable_testing()