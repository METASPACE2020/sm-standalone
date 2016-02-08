cmake_minimum_required(VERSION 2.8)
project(SM_standalone)

find_package(Boost REQUIRED COMPONENTS program_options)
include_directories(${Boost_INCLUDE_DIRS})

include(ExternalProject)

# MessagePack - simple serialization library
externalproject_add(msgpack
  GIT_REPOSITORY https://github.com/msgpack/msgpack-c
  GIT_TAG cpp-1.4.0
  CONFIGURE_COMMAND ""
  BUILD_COMMAND ""
  INSTALL_COMMAND ""
  ) # for C++ compilation is not required
ExternalProject_Get_Property(msgpack source_dir)
include_directories(${source_dir}/include)

# blosc - fast compressor for numeric data
externalproject_add(blosc
  GIT_REPOSITORY https://github.com/Blosc/c-blosc
  GIT_TAG v1.7.1
  INSTALL_COMMAND ""
  CMAKE_ARGS -DBUILD_TESTS:BOOL=OFF -DBUILD_BENCHMARKS:BOOL=OFF
  )

ExternalProject_Get_Property(blosc source_dir)
include_directories(${source_dir}/blosc/)

find_package(OpenMP REQUIRED)
find_package(LibXml2 REQUIRED)
find_package(ZLIB REQUIRED)
set(CMAKE_THREAD_PREFER_PTHREAD TRUE)
find_package(Threads REQUIRED)
include_directories(${LIBXML2_INCLUDE_DIR})

ExternalProject_Get_Property(blosc binary_dir)
set(BLOSC_LIBRARY ${binary_dir}/blosc/libblosc.a ${ZLIB_LIBRARY})

find_library(FFTW_LIBRARY NAMES fftw3 fftw)

if("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")
  message(FATAL_ERROR "In-source builds are strongly discouraged. "
    "Please remove CMakeCache.txt and the CMakeFiles/ directory, then build out-of-source.")
endif("${PROJECT_SOURCE_DIR}" STREQUAL "${PROJECT_BINARY_DIR}")

set(CMAKE_INCLUDE_CURRENT_DIR ON)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS} -g -std=c++11")
add_compile_options(-O2)

set(SRC_COMMON
  ims/image_measures.cpp
  ms/isocalc.cpp
  ms/isotope_pattern.cpp
  ms/periodic_table.cpp
  )

set(LIB_COMMON
  ${LIBXML2_LIBRARIES}
  ${FFTW_LIBRARY}
  ${BLOSC_LIBRARY}
  ${CMAKE_THREAD_LIBS_INIT}
  ${Boost_PROGRAM_OPTIONS_LIBRARY}
  )

add_executable(convert tools/convert.cpp
  imzml/reader.cpp ${SRC_COMMON})
target_link_libraries(convert ${LIB_COMMON})
add_dependencies(convert blosc)

add_executable(isocalc tools/isocalc.cpp ${SRC_COMMON})
target_link_libraries(isocalc ${LIB_COMMON})
add_dependencies(isocalc blosc msgpack)

add_library(imzb SHARED cffi/imzb_cffi.cpp)
target_link_libraries(convert ${LIBXML2_LIBRARIES} ${BLOSC_LIBRARY})
add_dependencies(imzb blosc)

set(CMAKE_EXPORT_COMPILE_COMMANDS 1) # for Vim autocompletion
