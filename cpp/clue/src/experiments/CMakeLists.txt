if(NOT TARGET ${PROJECT_NAME}-experiments)
  file(GLOB_RECURSE EXPERIMENTS_HEADERS ${PROJECT_SOURCE_DIR}/include/experiments/*.hpp)

  # add DD Package library
  add_library(
    ${PROJECT_NAME}-experiments
    ${EXPERIMENTS_HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/Experiment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/SATExperiment.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GroverExperiment.cpp
    )
  target_link_libraries(${PROJECT_NAME}-experiments PUBLIC ${PROJECT_NAME})
  add_library(${PROJECT_NAME}::Experiments ALIAS ${PROJECT_NAME}-experiments)
endif()
