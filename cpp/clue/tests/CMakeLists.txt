# Macro for doing the process for a new test file
macro(ADD_TEST_EXECUTABLE testname)
    add_executable(${PROJECT_NAME}_${testname}
                    ${CMAKE_CURRENT_SOURCE_DIR}/${testname}.cpp)
    target_link_libraries(${PROJECT_NAME}_${testname} PRIVATE ${PROJECT_NAME})
    add_test(NAME ${testname} COMMAND ${PROJECT_NAME}_${testname})
endmacro()

ADD_TEST_EXECUTABLE(LA_Subspace_Tests)
ADD_TEST_EXECUTABLE(LA_Vector_Tests)
ADD_TEST_EXECUTABLE(RF_Tests)