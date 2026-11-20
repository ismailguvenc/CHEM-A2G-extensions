# scm_mimo: 3GPP TR 25.996 SCM MIMO extension for CHEM (Nov 2026).
# Added by the top-level CMakeLists.txt when CHEM_EXT_SCM_MIMO=ON.
# Needs the a2g_doppler model library (geometry, trajectories, Doppler
# options), so CHEM_EXT_A2G_DOPPLER must also be ON.
# Targets:
#   chem_scm_model       static library, pure C++17 SCM model
#   scm_model_tests      unit tests of the model (GoogleTest)
#   scm_harness          offline MIMO harness: N IQ streams -> CHEM -> N IQ streams
#   chem_taps_timeline   writes a CIR timeline (a2g or scm model) for the
#                        live-twin driver (CHEM API, set_cir per port pair)
cmake_minimum_required(VERSION 3.15)

if(NOT TARGET chem_a2g_doppler_model)
  message(FATAL_ERROR "scm_mimo needs the a2g_doppler extension (CHEM_EXT_A2G_DOPPLER=ON)")
endif()

option(CHEM_EXT_SCM_MIMO_TESTS "Build scm_mimo unit tests" ON)
option(CHEM_EXT_SCM_MIMO_TOOLS "Build scm_mimo harness and timeline tool" ON)

add_library(chem_scm_model STATIC model/src/scm.cpp)
target_include_directories(chem_scm_model PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/model/include)
target_link_libraries(chem_scm_model PUBLIC chem_a2g_doppler_model)
target_compile_features(chem_scm_model PUBLIC cxx_std_17)
set_target_properties(chem_scm_model PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_compile_options(chem_scm_model PRIVATE -fno-fast-math -fcx-limited-range -O3)

if(TARGET CHEM)
  target_link_libraries(CHEM PRIVATE chem_scm_model)
  target_compile_definitions(CHEM PRIVATE CHEM_EXT_SCM_MIMO)
endif()
if(TARGET CHEM_LIB)
  target_link_libraries(CHEM_LIB PUBLIC chem_scm_model)
  target_compile_definitions(CHEM_LIB PUBLIC CHEM_EXT_SCM_MIMO)
endif()
if(TARGET achem_benchmark)
  target_link_libraries(achem_benchmark PRIVATE chem_scm_model)
  target_compile_definitions(achem_benchmark PRIVATE CHEM_EXT_SCM_MIMO)
endif()

if(CHEM_EXT_SCM_MIMO_TESTS)
  find_package(GTest)
  if(GTest_FOUND)
    enable_testing()
    add_executable(scm_model_tests tests/test_scm.cpp)
    target_link_libraries(scm_model_tests PRIVATE chem_scm_model GTest::gtest GTest::gtest_main)
    target_compile_options(scm_model_tests PRIVATE -fno-fast-math)
    add_test(NAME scm_model COMMAND scm_model_tests)
  endif()
endif()

if(CHEM_EXT_SCM_MIMO_TOOLS)
  if(TARGET CHEM_LIB)
    # uses the JSON parsers of both extensions (in CHEM_LIB)
    add_executable(chem_taps_timeline tools/chem_taps_timeline.cpp)
    target_link_libraries(chem_taps_timeline PRIVATE CHEM_LIB)
    target_include_directories(chem_taps_timeline PRIVATE
      ${PROJECT_SOURCE_DIR}/include ${chem_extension_include_dirs} ${PROJECT_SOURCE_DIR}/src)
    set_target_properties(chem_taps_timeline PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}")
    add_executable(chem_sionna_probe tools/chem_sionna_probe.cpp)
    target_link_libraries(chem_sionna_probe PRIVATE CHEM_LIB)
    target_include_directories(chem_sionna_probe PRIVATE
      ${PROJECT_SOURCE_DIR}/include ${chem_extension_include_dirs} ${PROJECT_SOURCE_DIR}/src)
    set_target_properties(chem_sionna_probe PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}")
    add_executable(scm_harness tools/scm_harness.cpp)
    target_link_libraries(scm_harness PRIVATE CHEM_LIB)
    target_include_directories(scm_harness PRIVATE
      ${PROJECT_SOURCE_DIR}/include ${chem_extension_include_dirs} ${PROJECT_SOURCE_DIR}/src)
    set_target_properties(scm_harness PROPERTIES
      RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin/${CMAKE_BUILD_TYPE}")
  endif()
endif()
