cmake_minimum_required(VERSION 3.0)
project(mc-pricer CXX)

# Prefix
if (NOT PREFIX)
    set(PREFIX ${CMAKE_CURRENT_BINARY_DIR})
endif (NOT PREFIX)
set(CMAKE_INSTALL_PREFIX ${PREFIX} CACHE STRING "Installation prefix." FORCE)

if (CXX)
    set(CMAKE_CXX_COMPILER ${CXX})
endif (CXX)

set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# Release or Debug
if (NOT CMAKE_BUILD_TYPE)
    message(STATUS "Setting build type to 'Debug' as none was specified.")
    set(CMAKE_BUILD_TYPE Debug CACHE STRING "Choose the type of build." FORCE)
endif ()


# Detect PNL
find_package(Pnl REQUIRED)
set(LIBS ${LIBS} ${PNL_LIBRARIES})
include_directories(${PNL_INCLUDE_DIRS})
if(${CMAKE_BUILD_TYPE} STREQUAL "Release")
    add_definitions(-DPNL_RANGE_CHECK_OFF)
endif()

# Add jlparser
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
if (WIN32)
    include(CheckCXXCompilerFlag)
    CHECK_CXX_COMPILER_FLAG("/std:c++latest" _cpp_latest_flag_supported)
    if (_cpp_latest_flag_supported)
        add_compile_options("/std:c++latest")
    endif()
endif()
add_subdirectory(src/3rdparty/jlparser)

### Ne rien modifier dans la partie au dessus de cette ligne.


# Test the parser
add_executable(test-parser src/test_parser.cpp)
target_link_libraries(test-parser jl_parser ${LIBS})

add_executable(mc-example src/mc_example.cpp)
target_link_libraries(mc-example jl_parser ${LIBS})

add_executable(Price0 src/Price0.cpp)
target_link_libraries(Price0 ${LIBS})

# add_executable(Hedge src/Hedge.cpp)
# target_link_libraries(Hedge ${LIBS})

###Google Tests
add_subdirectory(googletest)


enable_testing()

add_executable(
  hello_test
  src/hello_test.cpp
)
target_link_libraries(
  hello_test
  gtest_main
  ${LIBS}
)

add_executable(
  test_basket
  src/test_basket.cpp
)
target_link_libraries(
  test_basket
  gtest_main
  ${LIBS}
)

add_executable(
  test_asian
  src/test_asian.cpp
)
target_link_libraries(
  test_asian
  gtest_main
  ${LIBS}
)

add_executable(
  test_performance
  src/test_performance.cpp
)
target_link_libraries(
  test_performance
  gtest_main
  ${LIBS}
)

add_executable(
  test_portfolio
  src/test_portfolio.cpp
)
target_link_libraries(
  test_portfolio
  gtest_main
  ${LIBS}
)

include(GoogleTest)
gtest_discover_tests(hello_test)
pnl_add_postbuild(hello_test)
###


# Print compiler flags
get_directory_property(DirDefs COMPILE_DEFINITIONS)
message(STATUS "COMPILE_DEFINITIONS = ${DirDefs}")
if (${CMAKE_BUILD_TYPE} STREQUAL "Release")
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_RELEASE} ${CMAKE_LINKER_EXE_FLAGS}" )
else()
    message(STATUS "C++ Flags: ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_DEBUG} ${CMAKE_LINKER_EXE_FLAGS}" )
endif()
