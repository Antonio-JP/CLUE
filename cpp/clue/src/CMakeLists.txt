# add MQT::Core target
set(BUILD_MQT_DDSIM_TESTS
    OFF
    CACHE BOOL "Build MQT DDSIM tests")
add_subdirectory("${PROJECT_SOURCE_DIR}/extern/mqt-ddsim" "extern/mqt-ddsim"
                 EXCLUDE_FROM_ALL)


file(GLOB_RECURSE CLUE_HEADERS ${PROJECT_SOURCE_DIR}/include/*.hpp)

add_library(
    ${PROJECT_NAME}
    ${CLUE_HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/RationalFunctions.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Linalg_Vector.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Linalg_Subspace.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/Linalg_DD.cpp)

# set include directories
target_include_directories(${PROJECT_NAME} PUBLIC ${PROJECT_SOURCE_DIR}/include
                                                  ${PROJECT_BINARY_DIR}/include)


target_link_libraries(${PROJECT_NAME} PUBLIC MQT::ddsim)

# add alias
add_library(${PROJECT_NAME}::CPP ALIAS ${PROJECT_NAME})


# add experiments package library
add_subdirectory(experiments)