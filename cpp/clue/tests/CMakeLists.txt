project(CLUE_Tests LANGUAGES CXX)

add_executable(CLUE_RF_Tests src/RationalFunctions_Tests.cpp ../src/RationalFunctions.cpp)

# target_link_libraries(CLUE_RF_Tests PRIVATE MyLib)

add_test(NAME CLUE_RF COMMAND CLUE_RF_Tests)