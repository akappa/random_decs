project(matcher CXX)
set(PROJECT_VENDOR "Andrea Farruggia")
set(PROJECT_CONTACT "a.farruggia@di.unipi.it")
set(PROJECT_DESCRIPTION "Stand-alone Brotli matching strategy")

cmake_minimum_required(VERSION 2.8.7)
set(CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/CMakeModules")
include(Cpp14)
include(NinjaColorFix)
include(ReleaseDefault)
include(submodules)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_FASTDEB "-g -Wall -O3")
set(CMAKE_CXX_FLAGS_DEBUG "-g -Wall -O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -ffast-math -funroll-loops -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE_LTO "-O4 -march=native -ffast-math -funroll-loops -DNDEBUG")
set(CMAKE_CXX_FLAGS_TEST "-Wall -O3 -march=native -ffast-math -funroll-loops -fsanitize=address -DNDEBUG")

include_directories("${PROJECT_SOURCE_DIR}/include")

# Boost
find_package(Boost 1.46 REQUIRED program_options)
include_directories(${Boost_INCLUDE_DIR})
link_directories(${Boost_LIBRARY_DIRS})
list(APPEND AK_M_I_INCLUDE_DIRS "${Boost_INCLUDE_DIR}")

# Libraries
file(GLOB libraries libs/*.c*)
add_library(matcher_libs ${libraries})

add_subdirectory(external/match_serialize)
include_directories(${MatchSerializer_INCLUDE_DIRS})


# Main executables
function(exec_add binary)
  add_executable(${binary} ${binary}.cpp)
  target_link_libraries(${binary} matcher_libs match_serialize_libs ${Boost_LIBRARIES} ${ARGN})
endfunction()

exec_add(matcher)
