# a2g_doppler: per-tap A2G Doppler extension for CHEM.
# Added by the top-level CMakeLists.txt when CHEM_EXT_A2G_DOPPLER=ON.
# Targets:
#   chem_a2g_doppler_model   static library, pure C++17 model (no CHEM deps)
#   a2g_model_tests          unit tests of the model (GoogleTest)
#   a2g_hook_tests           tests through chem::Channel and the extension
#                            (needs -DCHEM_ENABLE_UNIT_TESTING=ON)
#   a2g_harness              offline harness: IQ file -> CHEM channel -> IQ file
cmake_minimum_required(VERSION 3.15)

option(CHEM_EXT_A2G_DOPPLER_TESTS "Build a2g_doppler unit tests" ON)
option(CHEM_EXT_A2G_DOPPLER_HARNESS "Build the a2g_doppler offline harness" ON)

add_library(chem_a2g_doppler_model STATIC model/src/model.cpp)
target_include_directories(chem_a2g_doppler_model PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/model/include)
target_compile_features(chem_a2g_doppler_model PUBLIC cxx_std_17)
set_target_properties(chem_a2g_doppler_model PROPERTIES POSITION_INDEPENDENT_CODE ON)
# CHEM builds with -ffast-math; the model keeps IEEE semantics. Without
# fast-math, GCC calls __muldc3 for every complex multiply; -fcx-limited-range
# inlines complex multiply (no inf/NaN recovery inside complex products).
target_compile_options(chem_a2g_doppler_model PRIVATE -fno-fast-math -fcx-limited-range -O3)

# Link the model into the CHEM targets and enable the extension source
if(TARGET CHEM)
  target_link_libraries(CHEM PRIVATE chem_a2g_doppler_model)
  target_compile_definitions(CHEM PRIVATE CHEM_EXT_A2G_DOPPLER)
endif()
if(TARGET CHEM_LIB)
  target_link_libraries(CHEM_LIB PUBLIC chem_a2g_doppler_model)
  target_compile_definitions(CHEM_LIB PUBLIC CHEM_EXT_A2G_DOPPLER)
endif()
if(TARGET achem_benchmark)
  target_link_libraries(achem_benchmark PRIVATE chem_a2g_doppler_model)
  target_compile_definitions(achem_benchmark PRIVATE CHEM_EXT_A2G_DOPPLER)
endif()

if(CHEM_EXT_A2G_DOPPLER_TESTS)
  find_package(GTest)
  if(GTest_FOUND)
    enable_testing()
    add_executable(a2g_model_tests tests/test_model.cpp)
    target_link_libraries(a2g_model_tests PRIVATE chem_a2g_doppler_model GTest::gtest GTest::gtest_main)
    target_compile_options(a2g_model_tests PRIVATE -fno-fast-math)
    add_test(NAME a2g_model COMMAND a2g_model_tests)
    if(TARGET CHEM_LIB)
      add_executable(a2g_hook_tests tests/test_hook.cpp)
      target_link_libraries(a2g_hook_tests PRIVATE CHEM_LIB GTest::gtest GTest::gtest_main)
      target_include_directories(a2g_hook_tests PRIVATE
        ${PROJECT_SOURCE_DIR}/include ${chem_extension_include_dirs} ${PROJECT_SOURCE_DIR}/src)
      add_test(NAME a2g_hook COMMAND a2g_hook_tests)
    endif()
  else()
    message(STATUS "a2g_doppler: GTest not found, unit tests skipped")
  endif()
endif()

if(CHEM_EXT_A2G_DOPPLER_HARNESS)
  if(TARGET CHEM_LIB)
    add_executable(a2g_harness tools/a2g_harness.cpp)
    target_link_libraries(a2g_harness PRIVATE CHEM_LIB)
  else()
    add_executable(a2g_harness tools/a2g_harness.cpp ${sources})
    target_link_libraries(a2g_harness PRIVATE ${privatelibs} chem_a2g_doppler_model)
    target_compile_definitions(a2g_harness PRIVATE CHEM_EXT_A2G_DOPPLER)
  endif()
  target_include_directories(a2g_harness PRIVATE
    ${PROJECT_SOURCE_DIR}/include ${chem_extension_include_dirs} ${PROJECT_SOURCE_DIR}/src)
  set_target_properties(a2g_harness PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}")
endif()
