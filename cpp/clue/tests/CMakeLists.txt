project(CLUE_Tests LANGUAGES CXX)

add_executable(CLUE_RF_Tests src/RationalFunctions_Tests.cpp ../src/RationalFunctions.cpp)
add_executable(CLUE_LA_Vector_Tests src/Linalg_Vector_Tests.cpp ../src/Linalg_Vector.cpp ../include/Linalg.hpp)
add_executable(CLUE_LA_Subspace_Tests src/Linalg_Subspace_Tests.cpp ../src/Linalg_Vector.cpp ../src/Linalg_Subspace.cpp ../include/Linalg.hpp)

# target_link_libraries(CLUE_RF_Tests PRIVATE MyLib)

add_test(NAME CLUE_RF COMMAND CLUE_RF_Tests)
add_test(NAME CLUE_LA_Vector COMMAND CLUE_LA_Vector_Tests)
add_test(NAME CLUE_LA_Subspace COMMAND CLUE_LA_Subspace_Tests)