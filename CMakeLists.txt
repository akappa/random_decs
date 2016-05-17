project(random_decs)
cmake_minimum_required(VERSION 2.8)

set(PROJECT_VENDOR "Andrea Farruggia")
set(PROJECT_CONTACT "a.farruggia@di.unipi.it")
set(PROJECT_DESCRIPTION "In-memory compressors/decompressors")

set(CMAKE_SOURCE_DIR ${PROJECT_SOURCE_DIR}/src)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${PROJECT_SOURCE_DIR}/cmake/Modules/")
set(LIB_DIR ${PROJECT_SOURCE_DIR}/libs)
set(LZ4_DIR ${PROJECT_SOURCE_DIR}/lz4)
set(CMAKE_PREFIX_PATH "~/.local/")

# C++11 support, Release by default
include(Cpp11)
include(NinjaColorFix)
include(ReleaseDefault)
include(submodules)
message(STATUS "Build type: ${CMAKE_BUILD_TYPE}")

set(CMAKE_CXX_FLAGS_FASTDEB "${GCC_DEBUG_FLAGS} -g -Wall -O3")
set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -g -Wall -O0")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -DNDEBUG -march=native")

set(CMAKE_C_FLAGS_FASTDEB "${GCC_DEBUG_FLAGS}  -Wall -W -Wundef -Wcast-align -Wno-implicit-function-declaration -Os -march=native -std=c99")
set(CMAKE_C_FLAGS_DEBUG "${GCC_DEBUG_FLAGS}  -Wall -W -Wundef -Wcast-align -Wno-implicit-function-declaration -Os -march=native -std=c99")
set(CMAKE_C_FLAGS_RELEASE "${GCC_DEBUG_FLAGS}  -Wall -W -Wundef -Wcast-align -Wno-implicit-function-declaration -Os -march=native -std=c99")

#include Brotli
add_subdirectory(brotli)
include_directories(brotli/include)
set(BROTLI_LIB "libbrotli")

# Include directory
include_directories("${PROJECT_SOURCE_DIR}/include")

# Boost
find_package(Boost 1.50 COMPONENTS iostreams program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_package(snappy)
if (NOT snappy_FOUND)
	message(FATAL_ERROR "Snappy not found")
	return()
endif()
include_directories(${snappy_INCLUDE_DIRS})

# Libraries
add_library(io ${LIB_DIR}/io.cpp)

set_source_files_properties(lz4.c lz4hc.c PROPERTIES COMPILE_FLAGS "-Wall -W -Wundef -Wcast-align -Wno-implicit-function-declaration -Os -march=native -std=c99")
add_library(lz4lc OBJECT ${LZ4_DIR}/lz4.c)
add_library(lz4hc OBJECT ${LZ4_DIR}/lz4hc.c)
add_library(lz4 $<TARGET_OBJECTS:lz4lc> $<TARGET_OBJECTS:lz4hc>)

# Booster
add_subdirectory(external/booster)

# Executable
add_executable(simple-lz4 ${CMAKE_SOURCE_DIR}/simple-lz4.cpp)
target_link_libraries(simple-lz4 io lz4)

add_executable(simple-lzma ${CMAKE_SOURCE_DIR}/simple-lzma.cpp)
target_link_libraries(simple-lzma io lzma)

add_executable(simple-snappy ${CMAKE_SOURCE_DIR}/simple-snappy.cpp)
target_link_libraries(simple-snappy ${snappy_LIBRARIES} io)

add_executable(simple-zlib ${CMAKE_SOURCE_DIR}/simple-zlib.cpp)
target_link_libraries(simple-zlib ${Boost_LIBRARIES} io)

add_executable(simple-bzip2 ${CMAKE_SOURCE_DIR}/simple-bzip2.cpp)
target_link_libraries(simple-bzip2 ${Boost_LIBRARIES} io)

add_executable(simple-brotli ${CMAKE_SOURCE_DIR}/simple-brotli.cpp)
target_link_libraries(simple-brotli ${Boost_LIBRARIES} ${BROTLI_LIB} io)
