# Macro to add a new executable and link the appropriate libraries
macro(ADD_APP_EXECUTABLE runname)
    add_executable(${PROJECT_NAME}_${runname}
                   ${CMAKE_CURRENT_SOURCE_DIR}/${runname}.cpp)
    target_link_libraries(${PROJECT_NAME}_${runname} PRIVATE ${PROJECT_NAME})
    target_link_libraries(${PROJECT_NAME}_${runname} PUBLIC ${ARGN})
endmacro()

ADD_APP_EXECUTABLE(main)
ADD_APP_EXECUTABLE(main_DD)
ADD_APP_EXECUTABLE(main_QASM)