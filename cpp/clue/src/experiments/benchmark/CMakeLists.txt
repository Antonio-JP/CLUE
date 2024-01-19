if(NOT TARGET ${PROJECT_NAME}-experiments-benchmark)
  file(GLOB_RECURSE EXPERIMENTS_BENCHMARK_HEADERS ${PROJECT_SOURCE_DIR}/include/experiments/benchmark/*.hpp)

  # add DD Package library
  add_library(
    ${PROJECT_NAME}-experiments-benchmark
    ${EXPERIMENTS_BENCHMARK_HEADERS}
    ${CMAKE_CURRENT_SOURCE_DIR}/AmplitudeEstimation.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/DeutschJozsa.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GHZ.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/GraphState.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/HHL.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PortfolioQAOA.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PortfolioVQE.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PricingCall.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/PricingPut.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QPEExact.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QPEInexact.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QuantumFourierTransform.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QuantumNeuralNetwork.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/QuantumWalk.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/TravellingSalesman.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/VariationalQuantumEigensolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/WState.cpp
    )
    target_link_libraries(${PROJECT_NAME}-experiments-benchmark PUBLIC ${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME}-experiments-benchmark PUBLIC ${PROJECT_NAME}-experiments)
  add_library(${PROJECT_NAME}::Experiments::Benchmark ALIAS ${PROJECT_NAME}-experiments-benchmark)
endif()